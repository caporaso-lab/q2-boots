# helpers
_copy-env-file:
	cp environment-files/readthedocs.yml book/_static/environment.yml

_copy-data:
	@if [ -d "book/data/" ]; then \
		cp -r book/data/ book/_build/html/data/; \
		echo "Copied data directory."; \
	else \
		echo "book/data/ not found, skipping copy."; \
	fi

_build-html:
	cd book && jupyter book build --html

_build-fast-preview:
	cd book && Q2DOC_FASTMODE= jupyter book build --html

_build-preview:
	cd book && Q2DOC_PREVIEW= jupyter book build --html

# main targets

autodoc:
	q2doc autodoc --singlepage --plugin boots --output plugin-reference book


html: _copy-env-file _build-html _copy-data

fast-preview: _copy-env-file _build-fast-preview _copy-data

preview: _copy-env-file _build-preview _copy-data

serve:
	npx serve book/_build/html/ -p 4000

clean:
	rm -rf book/_build/html/