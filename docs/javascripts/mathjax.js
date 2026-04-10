window.MathJax = {
  tex: {
    packages: { "[+]": ["ams"] },
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\[", "\\]"], ["$$", "$$"]],
    processEscapes: true,
    processEnvironments: true,
    macros: {
      // Provide a stable fallback for docs that use \boldsymbol{...}.
      // This keeps bold Greek/math symbols rendering even if the runtime
      // doesn't expose \boldsymbol directly.
      boldsymbol: ["\\pmb{#1}", 1]
    }
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
