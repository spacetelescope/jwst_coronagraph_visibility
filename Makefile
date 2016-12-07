app: clean
	pyinstaller run_jwst_coronagraph_visibility.spec

clean:
	find . -iname "*.pyc" -or -name "__pycache__" -delete
	rm -rf dist build

test:
	py.test
