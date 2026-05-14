// Cheap-trick client-side password gate. Password stored as SHA-256 hash
// rather than plaintext — keeps casual reads from picking up the password
// directly from this file. NOT real security:
//   - hash is visible; anyone can brute-force a 4-digit space in <1s
//   - sessionStorage can be set directly via devtools
//   - just keeps the URL un-shareable to random crawlers / link clicks
// Session persists until browser/tab close.
(function () {
    if (sessionStorage.getItem('ancestryAuth') === 'ok') return;

    var EXPECTED = '2926a2731f4b312c08982cacf8061eb14bf65c1a87cc5d70e864e079c6220731';
    var pw = prompt('Enter password to view:');

    // Hide everything immediately while async hash resolves.
    document.write('<style id="auth-block">body{display:none!important}</style>');

    function denied() {
        if (document.body) renderDenied();
        else document.addEventListener('DOMContentLoaded', renderDenied);
    }
    function renderDenied() {
        document.body.innerHTML =
            '<div style="position:fixed;inset:0;display:flex;align-items:center;justify-content:center;background:#1a1a1a;color:#999;font-family:-apple-system,sans-serif;z-index:99999;">' +
            '<div style="text-align:center"><h1 style="font-weight:300;font-size:1.4em">Access denied</h1>' +
            '<p style="margin-top:14px;color:#666;font-size:0.9em">Refresh to retry.</p></div></div>';
        document.body.style.cssText = 'margin:0;padding:0;display:block!important';
    }
    function allow() {
        sessionStorage.setItem('ancestryAuth', 'ok');
        var s = document.getElementById('auth-block');
        if (s) s.remove();
    }

    if (!pw) { denied(); return; }

    var buf = new TextEncoder().encode(pw);
    crypto.subtle.digest('SHA-256', buf).then(function (hash) {
        var hex = Array.from(new Uint8Array(hash))
            .map(function (b) { return b.toString(16).padStart(2, '0'); })
            .join('');
        if (hex === EXPECTED) allow(); else denied();
    }).catch(function () { denied(); });
})();
