"""HTTP server for the mirror with a DOI-style redirect."""
import http.server, sys, os, functools
root, port = sys.argv[1], int(sys.argv[2])
class H(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/doi/'):
            self.send_response(302); self.send_header('Location', f'http://127.0.0.1:{port}/records/14897628'); self.end_headers(); return
        if self.path.startswith('/records/'):
            self.send_response(200); self.end_headers(); self.wfile.write(b'record page'); return
        return super().do_GET()
    def log_message(self, *a): pass
http.server.ThreadingHTTPServer(('127.0.0.1', port), functools.partial(H, directory=root)).serve_forever()
