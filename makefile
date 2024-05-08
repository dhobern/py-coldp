reversion:
	python version.py increment

doc : 
	python version.py update docs/source/conf.py '^release = .*' 'release = "{version}"'
	python version.py update README.md 'documentation for COLDP version: [0-9\.]+' 'documentation for COLDP version: {version}'
	cd docs; ./make.bat html

package :
	python version.py update pyproject.toml 'version = .*' 'version = "{version}"'
	py -m pip install --upgrade build
	py -m build

build: doc package

publish : reversion build
	py -m twine upload --repository pypi dist/*`cat VERSION`*.*

clean :
	cd docs; ./make.bat clean
	rm -fR dist/*

