#!/bin/bash

# This script is used to publish the website for multiple cases to a remote host. 
# It reads the list of plot locations from a file called "plot_locs_list.txt" and uses scp to copy the website files 
# from each plot location to the specified remote host and directory.

# Usage: ./publish_multi_case.sh <remote_host> <remote_dir>
# NOTES: 
#   - The "plot_locs_list.txt" file should contain the paths to the plot locations, one per line.
#   - The path on the remote server has to already exist on the server side and the user must have write permissions to that path.

#local_dir=
remote_host="$1"
remote_dir="$2"
#while IFS= read -r item; do
#    echo "Processing $item"
#    scp -r "$item"/website/* "${remote_host}:${remote_dir}/"
#done < plot_locs_list.txt

#copy_cmd="scp -r"
if [[ "$mode" == "remote" ]]; then
    #copy_cmd="scp -r"
    prefix="${remote_host}:"
else
    #copy_cmd="cp -r"
    prefix=""
fi

while IFS= read -r item; do
    echo $prefix
    dest="${prefix}${remote_dir}"
    echo $dest
    scp -r "$item"/website/* "$dest"
done < plot_locs_list.txt




#!/bin/bash

mode="$1"

if [[ "$mode" == "remote" ]]; then
    remote_host="$2"
    remote_dir="$3"
elif [[ "$mode" == "local" ]]; then
    local_dir="$2"
else
    echo "Usage:"
    echo "  $0 remote <host> <remote_dir>"
    echo "  $0 local <local_dir>"
    exit 1
fi

while IFS= read -r item; do
    if [[ "$mode" == "remote" ]]; then
        dest="${remote_host}:${remote_dir}/${item}/"
        echo "Copying to remote: $dest"
        scp file.nc "$dest"
    else
        dest="${local_dir}/${item}/"
        echo "Copying locally: $dest"
        cp file.nc "$dest"
    fi
done < plot_locs_list.txt
