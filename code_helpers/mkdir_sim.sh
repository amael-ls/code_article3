#!/bin/bash

#### Check number of arguments
if [[ "$#" -ne 1 ]]
then
	echo "You must provide the txt file containing the simulation parameters"
	exit
fi

#### Read file and create directories if necessary
input="$1"
while IFS='' read -r line
do
	key=${line%=*}
	key=$(sed "s/ //g" <<< $key)
	if [[ "$key" =~ ^(species_path|initPath|summaryFilePath|popDynFilePath)$ ]]
	then
		value=${line##*=}
		value=$(sed "s/ //g" <<< $value)

		if [[ -d "$value" ]]
		then
			echo "$value already exists on your filesystem. Nothing was done"
		else
			echo "Creating directory $value"
			mkdir -p $value
		fi
	fi
done < "$input"
