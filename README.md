# q2-boots

q2-boots is a [QIIME 2](https://qiime2.org) plugin for microbiome rarefaction analysis that is currently in **beta** testing.

## Installation

### Install Prerequisites

[Miniconda](https://conda.io/miniconda.html) provides the `conda` environment and package manager, and is currently the only supported way to install QIIME 2.
Follow the instructions for downloading and installing Miniconda.

After installing Miniconda and opening a new terminal, make sure you're running the latest version of `conda`, and install `wget` in case you don't already have it.

```bash
conda update conda
conda install wget
```

### QIIME 2024.10 development version of q2-boots

Installing the most recent development version of q2-boots allows you to access the most recent functionality, including some that depends on features being introduced in QIIME 2 2024.10.
It is hard to unambiguously reference development versions of software in publications however, so you'll likely want to re-run your boots analyses with a release version prior to publication.

```shell
conda env create -n q2-boots-2024.10 -f https://raw.githubusercontent.com/qiime2/q2-boots/main/environments/q2-boots-qiime2-amplicon-2024.10.yml
```

```shell
conda activate q2-boots-2024.10
```

## Testing and using `q2-boots`

After completing the install steps above, confirm that everything is working as expected by running:

```shell
make test
```

You should get a report that tests were run, and you should see that all tests passed and none failed.
It's usually ok if some warnings are reported.

If all of the tests pass, you're ready to use the plugin.
Start by making QIIME 2's command line interface aware of `q2-boots` by running:

```shell
qiime dev refresh-cache
```

You should then see the plugin in the list of available plugins if you run:

```shell
qiime info
```

You should be able to review the help text by running:

```shell
qiime boots --help
```

Have fun! 😎

## Getting help

If you need help with q2-boots, please get in touch on the [QIIME 2 Forum](https://forum.qiime2.org).

## About

The `q2-boots` Python package was built using content from the [QIIME 2 plugin template](https://develop.qiime2.org/en/latest/plugins/tutorials/create-from-template.html).
To learn more about `q2-boots`, refer to the [project website](https://github.com/qiime2/q2-boots).
To learn how to use QIIME 2, refer to the [QIIME 2 User Documentation](https://docs.qiime2.org).
To learn QIIME 2 plugin development, refer to [*Developing with QIIME 2*](https://develop.qiime2.org).

`q2-boots` is a QIIME 2 community plugin, meaning that it is not necessarily developed and maintained by the developers of QIIME 2.
Please be aware that because community plugins are developed by the QIIME 2 developer community, and not necessarily the QIIME 2 developers themselves, some may not be actively maintained or compatible with current release versions of the QIIME 2 distributions.
More information on development and support for community plugins can be found [here](https://library.qiime2.org).

![Boots covered in microorganisms](./images/q2-boots-ai-art.png)