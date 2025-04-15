# q2-boots Documentation source

[![Copier](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/copier-org/copier/master/img/badge/badge-grayscale-inverted-border-orange.json)](https://github.com/copier-org/copier)

## Development instructions

The following sub-sections illustrate how to develop this documentation.

### Create the development environment

To build this documentation locally for development purposes, first create your development environment.

```
cd q2-boots-docs
conda env create -n q2-boots-docs --file environment-files/readthedocs.yml
conda activate q2-boots-docs
q2doc refresh-cache
```


### Autogenerate plugin documentation

Then, create the auto-generated documentation.
This will only need to be run the first time you build your documentation, or something has changed in the plugin.

```
make autodoc
```


### Build the book

Next, build the book:

```
make html
```

(Alternatively, `make preview` or `make fast-preview` can speed up test builds.)

### Serve the book locally

Finally, run the following to serve the built documentation locally:

```
make serve
```

Have fun!
