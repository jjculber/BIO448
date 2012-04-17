#!/bin/bash
make clean
make
jar cvfm CGContent.jar manifest.txt bioGUI/*
