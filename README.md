# Pi0_analysis
Pi0 analysis code for reconstructing the Pi0 4-vectors from NPS

## Versioning

- Update `VERSION` for each release tag (example: `0.1.1`).
- Tag in git with the same value (example: `git tag pi0-pipeline-v0.1.1 && git push origin pi0-pipeline-v0.1.1`).

## Docker image build

```bash
export PI0_VERSION=$(cat VERSION)
docker build \
	--build-arg PI0_ANALYSIS_VERSION="$PI0_VERSION" \
	--build-arg VCS_REF="$(git rev-parse --short HEAD)" \
	--build-arg BUILD_DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
	-t your-registry/pi0-analysis:"$PI0_VERSION" .
```
