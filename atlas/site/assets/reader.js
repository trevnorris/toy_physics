
const filter = document.querySelector('[data-topic-filter]');
if (filter) {
  filter.addEventListener('input', () => {
    const query = filter.value.trim().toLowerCase();
    document.querySelectorAll('[data-topic-card]').forEach((card) => {
      card.hidden = query && !card.textContent.toLowerCase().includes(query);
    });
  });
}
