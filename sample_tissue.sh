#!/bin/bash

# This script is used to extract the sample type from the XML for the Castillo Fernandes data

cat `find . | grep sample.xml`  \
| egrep '(blood|namespace)' \
| paste -s \
| sed 's/>/\n/g' \
| egrep '(^W|^C|^wh|^co)' \
| sed 's/</\t/' \
| cut -f1 \
| paste - - \
| sort -u

