#!/bin/bash

#path to executable file that handles audio files
executable_path="/home/qwerty/FindSilence/findsilence"

#setting default values to parameters
sort_key="-k 1,1 -k 2,2"
threshold=-1
duration=0.25

args=( "$@" )

isfloat ()
{
	regexp='^[-+]?[0-9]*\.?[0-9]+$'

	if [[ $1 =~ $regexp ]]
	then
		return 1
	else
		return 0
	fi
}

for (( i=0; i < $#; ++i ))
do
	case ${args[i]:0:2} in
	"-t")
		isfloat ${args[++i]}

		if [ $? == 0 ]
		then
			echo "Wrong threshold value -- ${args[i]}"
			exit -1
		fi

		threshold=${args[i]}
		;;
	"-m")
		isfloat ${args[++i]}

		if [ $? == 0 ]
		then
			echo "Wrong minimun silence duration value -- ${args[i]}"
			exit -1
		fi

		duration=${args[i]}
		;;

	"-s")
		case ${args[i]:2:3} in
		'c')
			sort_key="-k 1,1 -k 2,2"
			;;
		't')
			sort_key="-k 2,2"
			;;
		*)
			echo "Wrong sorting key -- ${args[i]:2:3}" >&2
			exit -1
			;;
		esac
		;;
	*)
		echo "Wrong argument -- ${args[i]:0:2}" >&2
		exit -1
		;;
	esac
done

while read fname
do
	sox $fname -e signed-integer /tmp/fsilence_$$.wav
	echo $fname
#	$executable_path "/tmp/fsilence_$$.wav" $threshold $duration | sort -n $sort_key
	$executable_path "/tmp/fsilence_$$.wav" $threshold $duration
done

rm "/tmp/fsilence_$$.wav"
