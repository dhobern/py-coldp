reversion:
	python version.py increment

doc : 
	python version.py update docs/source/conf.py '^release = .*' 'release = "{version}"'
	python version.py update README.md 'documentation for COLDP version: [0-9\.]+' 'documentation for COLDP version: {version}'
	cd docs; make html; make simplepdf

package :
	python version.py update pyproject.toml 'version = .*' 'version = "{version}"'
	python -m pip install --upgrade build
	python -m build

build: doc package

publish : reversion clean build
	python -m twine upload --repository pypi dist/*`cat VERSION`*.*

clean :
	cd docs; make clean
	rm -fR dist/*

