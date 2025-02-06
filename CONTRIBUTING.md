# Contributing

If you would like to contribute to the development of `biosynfoni`, please open an issue or a pull request. We welcome contributions from the community.

## Development installation

Run the following command from the root of the project to install the project for development:

```bash
pip install -e .[dev]
```

Tests will be run automatically when you push a pull request. You can also run the tests locally with the following command:

```bash
pytest tests/
```

Please use `black` to format your code before submitting a pull request:

```bash
black .
```
