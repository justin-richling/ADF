#!/bin/bash

git fetch origin $1
git branch $1 FETCH_HEAD
git checkout $1
