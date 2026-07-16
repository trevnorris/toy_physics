#define _GNU_SOURCE
#include <dlfcn.h>
#include <elf.h>
#include <errno.h>
#include <fcntl.h>
#include <link.h>
#include <linux/openat2.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>

/* A self-contained SHA-256 keeps the auditor from adding a crypto-library
 * dependency before the loader audit callbacks become live. */
typedef struct { uint32_t h[8]; uint64_t bits; unsigned char block[64]; size_t used; } sha256_ctx;
static const uint32_t K[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};
static uint32_t rotr(uint32_t x, unsigned n) { return (x >> n) | (x << (32-n)); }
static void transform(sha256_ctx *c, const unsigned char *p) {
  uint32_t w[64],a,b,d,e,f,g,h,t1,t2,cc; unsigned i;
  for(i=0;i<16;i++) w[i]=((uint32_t)p[4*i]<<24)|((uint32_t)p[4*i+1]<<16)|((uint32_t)p[4*i+2]<<8)|p[4*i+3];
  for(i=16;i<64;i++){uint32_t s0=rotr(w[i-15],7)^rotr(w[i-15],18)^(w[i-15]>>3),s1=rotr(w[i-2],17)^rotr(w[i-2],19)^(w[i-2]>>10);w[i]=w[i-16]+s0+w[i-7]+s1;}
  a=c->h[0];b=c->h[1];cc=c->h[2];d=c->h[3];e=c->h[4];f=c->h[5];g=c->h[6];h=c->h[7];
  for(i=0;i<64;i++){uint32_t S1=rotr(e,6)^rotr(e,11)^rotr(e,25),ch=(e&f)^((~e)&g),S0=rotr(a,2)^rotr(a,13)^rotr(a,22),maj=(a&b)^(a&cc)^(b&cc);t1=h+S1+ch+K[i]+w[i];t2=S0+maj;h=g;g=f;f=e;e=d+t1;d=cc;cc=b;b=a;a=t1+t2;}
  c->h[0]+=a;c->h[1]+=b;c->h[2]+=cc;c->h[3]+=d;c->h[4]+=e;c->h[5]+=f;c->h[6]+=g;c->h[7]+=h;
}
static void sha_init(sha256_ctx *c){uint32_t init[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};memcpy(c->h,init,sizeof init);c->bits=0;c->used=0;}
static void sha_update(sha256_ctx *c,const unsigned char *p,size_t n){c->bits+=(uint64_t)n*8;while(n){size_t take=64-c->used;if(take>n)take=n;memcpy(c->block+c->used,p,take);c->used+=take;p+=take;n-=take;if(c->used==64){transform(c,c->block);c->used=0;}}}
static void sha_final(sha256_ctx *c,unsigned char out[32]){unsigned i;c->block[c->used++]=0x80;if(c->used>56){while(c->used<64)c->block[c->used++]=0;transform(c,c->block);c->used=0;}while(c->used<56)c->block[c->used++]=0;for(i=0;i<8;i++)c->block[63-i]=(unsigned char)(c->bits>>(8*i));transform(c,c->block);for(i=0;i<8;i++){out[4*i]=c->h[i]>>24;out[4*i+1]=c->h[i]>>16;out[4*i+2]=c->h[i]>>8;out[4*i+3]=c->h[i];}}

static __thread int guard;
static int logfd = -2;
typedef struct { dev_t device; ino_t inode; unsigned state; } seen_file;
#define SEEN_CAPACITY 262144U
#define OBJECT_CAPACITY 32768U
static seen_file seen_files[SEEN_CAPACITY];
static uint64_t seen_hash(dev_t device,ino_t inode){
  uint64_t x=(uint64_t)device^((uint64_t)inode+0x9e3779b97f4a7c15ULL+((uint64_t)device<<6)+((uint64_t)device>>2));
  x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);
}
static int seen_before(dev_t device,ino_t inode){
  unsigned probe;
  uint64_t start=seen_hash(device,inode)&(SEEN_CAPACITY-1U);
  for(probe=0;probe<SEEN_CAPACITY;probe++){
    seen_file *slot=&seen_files[(start+probe)&(SEEN_CAPACITY-1U)];
    unsigned state=__atomic_load_n(&slot->state,__ATOMIC_ACQUIRE);
    while(state==1U)state=__atomic_load_n(&slot->state,__ATOMIC_ACQUIRE);
    if(state==2U){if(slot->device==device&&slot->inode==inode)return 1;continue;}
    if(__sync_bool_compare_and_swap(&slot->state,0U,1U)){
      slot->device=device;slot->inode=inode;__atomic_store_n(&slot->state,2U,__ATOMIC_RELEASE);return 0;
    }
  }
  return 1;
}
static int raw_open(const char *path,int flags){return (int)syscall(SYS_openat,AT_FDCWD,path,flags,0600);}
static int get_logfd(void){const char *p;if(logfd!=-2)return logfd;p=getenv("B2_FIRST_USE_LOG");if(!p||!*p){logfd=-1;return -1;}logfd=raw_open(p,O_WRONLY|O_CREAT|O_APPEND|O_CLOEXEC);return logfd;}
static int same_log_path(const char *path){const char *p=getenv("B2_FIRST_USE_LOG");return p&&path&&strcmp(p,path)==0;}
static void record_fd(int fd,const char *kind,const char *fallback){
  struct stat st; sha256_ctx c; unsigned char buf[65536],sum[32]; char hex[65],path[4096],proc[64],line[4600]; ssize_t n; off_t off=0; int out,i;
  if(guard||fd<0){return;} guard=1;
  if(fstat(fd,&st)||!S_ISREG(st.st_mode)){guard=0;return;}
  if(seen_before(st.st_dev,st.st_ino)){guard=0;return;}
  snprintf(proc,sizeof proc,"/proc/self/fd/%d",fd);n=readlink(proc,path,sizeof(path)-1);if(n>0)path[n]=0;else snprintf(path,sizeof path,"%s",fallback?fallback:"<unknown>");
  if(same_log_path(path)){guard=0;return;}
  sha_init(&c);while((n=pread(fd,buf,sizeof buf,off))>0){sha_update(&c,buf,(size_t)n);off+=n;}if(n<0){guard=0;return;}sha_final(&c,sum);for(i=0;i<32;i++)sprintf(hex+2*i,"%02x",sum[i]);hex[64]=0;
  out=get_logfd();if(out>=0){int len=snprintf(line,sizeof line,"%ld\t%s\t%s\t%llu\t%llu\t%s\n",(long)getpid(),kind,hex,(unsigned long long)st.st_dev,(unsigned long long)st.st_ino,path);if(len>0)syscall(SYS_write,out,line,(size_t)len);}
  guard=0;
}
static void record_path(const char *path,const char *kind){int fd;if(!path||!*path)path="/proc/self/exe";fd=raw_open(path,O_RDONLY|O_CLOEXEC);if(fd>=0){record_fd(fd,kind,path);syscall(SYS_close,fd);}}

/* LD_PRELOAD-compatible inventory for executables that reject an rtld audit
 * namespace (the host-licensed Wolfram kernel is one).  The constructor runs
 * before application main, and each subsequent dlopen/dlmopen is inventoried
 * before the call returns to its consumer. */
static __thread int snapshot_guard;
static uintptr_t seen_objects[OBJECT_CAPACITY]; static unsigned seen_object_count;
static int snapshot_one(struct dl_phdr_info *info,size_t size,void *data){const char *p;uintptr_t key;unsigned count,j;(void)size;(void)data;p=info&&info->dlpi_name?info->dlpi_name:"";key=(uintptr_t)(info?info->dlpi_addr:0)^(uintptr_t)(p&&*p?1:0);count=seen_object_count;if(count>OBJECT_CAPACITY)count=OBJECT_CAPACITY;for(j=0;j<count;j++)if(seen_objects[j]==key)return 0;j=__sync_fetch_and_add(&seen_object_count,1);if(j>=OBJECT_CAPACITY)return 0;seen_objects[j]=key;record_path(p&&*p?p:"/proc/self/exe",p&&*p?"LOAD":"EXEC");return 0;}
static void snapshot_objects(void){if(snapshot_guard)return;snapshot_guard=1;dl_iterate_phdr(snapshot_one,NULL);snapshot_guard=0;}
__attribute__((constructor)) static void first_use_startup(void){snapshot_objects();unsetenv("LD_AUDIT");}

unsigned int la_version(unsigned int v){return v>LAV_CURRENT?LAV_CURRENT:v;}
unsigned int la_objopen(struct link_map *map,Lmid_t lmid,uintptr_t *cookie){(void)lmid;(void)cookie;record_path(map&&map->l_name&&*map->l_name?map->l_name:"/proc/self/exe",map&&map->l_name&&*map->l_name?"LOAD":"EXEC");return 0;}

typedef int (*open_fn)(const char*,int,...); typedef int (*openat_fn)(int,const char*,int,...);
static int finish_open(int fd,int flags,const char *path){if(fd>=0 && (flags&O_ACCMODE)!=O_WRONLY)record_fd(fd,"OPEN",path);return fd;}
int open(const char *p,int f,...){static open_fn real;mode_t m=0;if(f&O_CREAT){va_list a;va_start(a,f);m=(mode_t)va_arg(a,int);va_end(a);}if(!real)real=(open_fn)dlsym(RTLD_NEXT,"open");return finish_open((f&O_CREAT)?real(p,f,m):real(p,f),f,p);}
int open64(const char *p,int f,...){static open_fn real;mode_t m=0;if(f&O_CREAT){va_list a;va_start(a,f);m=(mode_t)va_arg(a,int);va_end(a);}if(!real)real=(open_fn)dlsym(RTLD_NEXT,"open64");return finish_open((f&O_CREAT)?real(p,f,m):real(p,f),f,p);}
int openat(int d,const char *p,int f,...){static openat_fn real;mode_t m=0;if(f&O_CREAT){va_list a;va_start(a,f);m=(mode_t)va_arg(a,int);va_end(a);}if(!real)real=(openat_fn)dlsym(RTLD_NEXT,"openat");return finish_open((f&O_CREAT)?real(d,p,f,m):real(d,p,f),f,p);}
int openat64(int d,const char *p,int f,...){static openat_fn real;mode_t m=0;if(f&O_CREAT){va_list a;va_start(a,f);m=(mode_t)va_arg(a,int);va_end(a);}if(!real)real=(openat_fn)dlsym(RTLD_NEXT,"openat64");return finish_open((f&O_CREAT)?real(d,p,f,m):real(d,p,f),f,p);}
int openat2(int d,const char *p,const struct open_how *how,size_t size){int fd=(int)syscall(SYS_openat2,d,p,how,size);return finish_open(fd,how?(int)how->flags:O_RDONLY,p);}

typedef void *(*dlopen_fn)(const char*,int); typedef void *(*dlmopen_fn)(Lmid_t,const char*,int);
void *dlopen(const char *p,int f){static dlopen_fn real;void *h;if(!real)real=(dlopen_fn)dlsym(RTLD_NEXT,"dlopen");h=real(p,f);if(h)snapshot_objects();return h;}
void *dlmopen(Lmid_t ns,const char *p,int f){static dlmopen_fn real;void *h;if(!real)real=(dlmopen_fn)dlsym(RTLD_NEXT,"dlmopen");h=real(ns,p,f);if(h)snapshot_objects();return h;}
