#!/bin/bash

#path to executable file that handles audio files
executable_path="/home/qwerty/FindSilence/findsilence"

#setting default values to parameters
sort_key="-k 1,1 -k 2,2"
threshold=-1
si_duration=0.25
so_duration=0.0625

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
	case ${args[i]} in
	"--help")
		printf "Usage: findsilence [OPTIONS]\n\n"
		printf "  -t FLOAT\tset threshold value (automaticly by default)\n" 
		printf "  -ms FLOAT\tset minimum silence duration value (default value is %f)\n" $si_duration
		printf "  -mS FLOAT\tset minimum sound duration value (default value is %f)\n" $so_duration
		printf "  -st\t\tsort by time, when event occurred\n"
		printf "  -sc\t\tsort by channel, where event occurred (default)\n\n"
		printf "Example:\n"
		printf "  $ findsilence -t 0.5 -m 0.25 -st\n\n"

		exit 0
		;;
	"-t")
		isfloat ${args[++i]}

		if [ $? == 0 ]
		then
			echo "Wrong threshold value -- ${args[i]}"
			exit 1
		fi

		threshold=${args[i]}
		;;
	"-ms"|"-mS")
	
		isfloat ${args[i+1]}
		
		if [ $? == 0 ]
			then
				echo "Wrong minimun duration value -- ${args[i+1]}"
				exit 1
			fi

		case ${args[i]:2:3} in
		's')
			si_duration=${args[i+1]}
			;;
		'S')
			so_duration=${args[i+1]}		
			;;
		*)
			echo "Wrong argument -- ${args[i]}" >&2
			;;
		esac

		((++i))
		;;
		

	"-sc"|"-st")
		case ${args[i]:2:3} in
		'c')
			sort_key="-k 1,1 -k 2,2"
			;;
		't')
			sort_key="-k 2,2"
			;;
		*)
			echo "Wrong argument -- ${args[i]}" >&2
			exit 1
			;;
		esac
		;;
	*)
		echo "Wrong argument -- ${args[i]}" >&2
		exit -1
		;;
	esac
done

while true
do
	read fname
	tmp=$?
	
	if [[ "$fname" = "" && ${tmp} != 0 ]]
	then
		break
	fi
		
	sox $fname -e signed-integer /tmp/fsilence_$$.wav
	echo "F $fname"
	$executable_path "/tmp/fsilence_$$.wav" $threshold $si_duration $so_duration | sort -n $sort_key
done

if [ -f "/tmp/fsilence_$$.wav" ]
then
	rm "/tmp/fsilence_$$.wav"
fi
