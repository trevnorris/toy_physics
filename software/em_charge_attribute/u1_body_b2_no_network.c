#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <netinet/in.h>
#include <stddef.h>
#include <string.h>
#include <sys/socket.h>

/* Traced-mode production policy: local AF_UNIX sockets remain available to the
 * interpreter/license machinery; internet families are denied before a file
 * descriptor can be created or connected.  The runner additionally audits the
 * network syscall trace, so direct-syscall bypasses are observable. */

typedef int (*socket_fn)(int, int, int);
typedef int (*connect_fn)(int, const struct sockaddr *, socklen_t);
typedef int (*bind_fn)(int, const struct sockaddr *, socklen_t);
typedef ssize_t (*sendto_fn)(int, const void *, size_t, int, const struct sockaddr *, socklen_t);

static int external_address(const struct sockaddr *address) {
    if (address == NULL) return 0;
    if (address->sa_family == AF_INET) {
        const struct sockaddr_in *v4 = (const struct sockaddr_in *)address;
        return (ntohl(v4->sin_addr.s_addr) >> 24) != 127;
    }
    if (address->sa_family == AF_INET6) {
        const struct sockaddr_in6 *v6 = (const struct sockaddr_in6 *)address;
        return !IN6_IS_ADDR_LOOPBACK(&v6->sin6_addr);
    }
    return 0;
}

int socket(int domain, int type, int protocol) {
    static socket_fn real_socket = NULL;
    if (real_socket == NULL) real_socket = (socket_fn)dlsym(RTLD_NEXT, "socket");
    if (real_socket == NULL) {
        errno = ENOSYS;
        return -1;
    }
    return real_socket(domain, type, protocol);
}

int connect(int fd, const struct sockaddr *address, socklen_t length) {
    static connect_fn real_connect = NULL;
    if (external_address(address)) {
        errno = EACCES;
        return -1;
    }
    if (real_connect == NULL) real_connect = (connect_fn)dlsym(RTLD_NEXT, "connect");
    if (real_connect == NULL) {
        errno = ENOSYS;
        return -1;
    }
    return real_connect(fd, address, length);
}

int bind(int fd, const struct sockaddr *address, socklen_t length) {
    static bind_fn real_bind = NULL;
    struct sockaddr_storage local;
    const struct sockaddr *effective = address;
    if (address != NULL && address->sa_family == AF_INET && length >= sizeof(struct sockaddr_in)) {
        memcpy(&local, address, sizeof(struct sockaddr_in));
        ((struct sockaddr_in *)&local)->sin_addr.s_addr = htonl(INADDR_LOOPBACK);
        effective = (const struct sockaddr *)&local;
    } else if (address != NULL && address->sa_family == AF_INET6 && length >= sizeof(struct sockaddr_in6)) {
        memcpy(&local, address, sizeof(struct sockaddr_in6));
        ((struct sockaddr_in6 *)&local)->sin6_addr = in6addr_loopback;
        effective = (const struct sockaddr *)&local;
    }
    if (real_bind == NULL) real_bind = (bind_fn)dlsym(RTLD_NEXT, "bind");
    if (real_bind == NULL) {
        errno = ENOSYS;
        return -1;
    }
    return real_bind(fd, effective, length);
}

ssize_t sendto(int fd, const void *buffer, size_t size, int flags, const struct sockaddr *destination, socklen_t length) {
    static sendto_fn real_sendto = NULL;
    if (external_address(destination)) {
        errno = EACCES;
        return -1;
    }
    if (real_sendto == NULL) real_sendto = (sendto_fn)dlsym(RTLD_NEXT, "sendto");
    if (real_sendto == NULL) {
        errno = ENOSYS;
        return -1;
    }
    return real_sendto(fd, buffer, size, flags, destination, length);
}
