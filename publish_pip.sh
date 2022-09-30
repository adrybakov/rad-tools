#! /bin/zsh

rm -r dist
rm -r build
rm -r rad_tools.egg-info
python3 setup.py sdist bdist_wheel 
python3 -m twine upload --repository pypi dist/*