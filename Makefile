.PHONY: all lint test install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

install: all
	pip install .

dev: all
	pip install -e .

clean: distclean

distclean: ;

# doc build targets and helpers
_copy-env-file:
	cp docs/environment-files/readthedocs.yml docs/book/_static/environment.yml

_copy-data:
	@if [ -d "docs/book/data/" ]; then \
		cp -r docs/book/data/ docs/book/_build/html/data/; \
		echo "Copied data directory."; \
	else \
		echo "docs/book/data/ not found, skipping copy."; \
	fi

_build-html:
	cd docs/book && jupyter book build --html

_build-fast-preview:
	cd docs/book && Q2DOC_FASTMODE= jupyter book build --html

_build-preview:
	cd docs/book && Q2DOC_PREVIEW= jupyter book build --html

autodoc:
	cd docs/book && q2doc autodoc --singlepage --plugin boots --output plugin-reference .
	ls -al

html: _copy-env-file _build-html _copy-data

fast-preview: _copy-env-file _build-fast-preview _copy-data

preview: _copy-env-file _build-preview _copy-data

serve:
	npx serve docs/book/_build/html/ -p 4000

clean:
	rm -rf docs/book/_build/html/