#!/usr/bin/env bash

function HELP {
  echo -e \\n"Builds a txtable for tu_selecter from at gtf file. Writes to stdout."
  echo -e "Basic usage: build_txtable.sh gene_annotations.gtf"\\n
  echo -e "\t-h Displays this help message. No further functions are performed."
  exit 1
}

while getopts "h" opt; do
  case $opt in
      h)  #show help
      HELP
      ;;
   \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

gtf=${@:$OPTIND:1}

awk 'BEGIN {OFS="\t"}NR > 5 {if($3=="transcript"){gsub(/\"|\;/,"",$10);gsub(/\"|\;/,"",$12);gsub(/\"|\;/,"",$14);gsub(/\"|\;/,"",$18);print $1,$4,$5,$10,$12,$7,$18,$14}}' $gtf