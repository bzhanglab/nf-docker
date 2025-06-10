# MS²Rescore with Percolator Docker Image

This Docker image provides [MS²Rescore](https://github.com/CompOmics/ms2rescore) with [Percolator](https://github.com/percolator/percolator) as the rescoring engine.

## Features

- **MS²Rescore v3.1+**: Latest version with full feature support
- **Percolator v3.06.05**: Pre-compiled and ready to use
- **All dependencies included**: MS²PIP, DeepLC, and other feature generators
- **Multi-format support**: Works with various search engine outputs

## Building the Image

```bash
cd ms2rescore
docker build -t registry.gitlab.com/bzhanglab/ms2rescore:3.1.5 .
```

## Usage

```bash
docker run --rm -v /path/to/your/data:/data \
  registry.gitlab.com/bzhanglab/ms2rescore:3.1.5 \
  --config_file /data/config.toml
```

## Testing

Run the comprehensive test suite:

```bash
chmod +x test-docker.sh
./test-docker.sh
```

## Configuration

Use `percolator-config.toml` as a template for your configuration file.

## Supported Formats

- MaxQuant (`msms.txt`)
- Sage (`.sage.tsv`)
- MSGFPlus (`.mzid`)
- X!Tandem (`.xml`)
- And more... 