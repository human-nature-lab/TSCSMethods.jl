// Load Mermaid.js from CDN and initialize
(function() {
    // Load Mermaid.js from CDN
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/mermaid@10.6.1/dist/mermaid.min.js';
    script.onload = function() {
        // Initialize Mermaid with configuration
        mermaid.initialize({
            startOnLoad: true,
            theme: 'default',
            themeVariables: {
                fontFamily: 'Arial, sans-serif',
                fontSize: '14px'
            }
        });
    };
    document.head.appendChild(script);
})();