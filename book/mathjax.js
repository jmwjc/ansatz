// MathJax 3 configuration: enable $...$ inline and $$...$$ display delimiters.
// MathJax 2 (mdBook's built-in) does not enable $ inline math by default,
// so we load MathJax 3 with explicit delimiters for a nicer CHTML output.
window.MathJax = {
  tex: {
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$'], ['\\[', '\\]']],
    processEscapes: true
  },
  options: {
    skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
  }
};

(function () {
  var s = document.createElement('script');
  s.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js';
  s.async = true;
  document.head.appendChild(s);
})();
