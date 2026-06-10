#!/bin/sh

git clone --branch 1.1.11 --depth 1 https://github.com/SergiusTheBest/plog.git plog    
cd plog
rm -rf .circleci .git samples
rm -rf .appveyor.yml .editorconfig .gitignore .travis.yml
