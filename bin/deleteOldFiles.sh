#!/bin/bash

# delete files older than 15 days
find $HOME/$USER/Downloads/* -mtime +15 -exec rm -rf {} \;
find /tmp/*  -mtime +15 -exec rm -rf {} \;
