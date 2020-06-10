## START SCRIPT     START SCRIPT     START SCRIPT     START SCRIPT     START SCRIPT
## A script to sort and extract meaningful information from GapMinder Data
## Here: to extract for the year 2007 the 10 countries with the highest GDP per capita and the 10 countries with the lowest GDP per capita  and to sort them according to GDP
## to begin grab all data for the year 2007 from gapminder.txt
grep "\b2007\b" /mnt/c/users/sreev/KEYS/KEYS_Training/Data/Data_Workfolder/gapminder_by_country/gapminder.txt > /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/2007_data.txt
## ensure that there will be no confusion about columns by replacing empty spaces with '_'
sed -i 's/ /_/g' /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/2007_data.txt
## sort the 2007 data for all 142 countries by GDP in ascendig order
sort -k6 -b -n /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/2007_data.txt > /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/2007_sortedbyGDP.txt
## extract the 10 richest countries from the end of the file and present them in reverse order, i.e., the country with the highest GDP on top
tail -10 /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/2007_sortedbyGDP.txt | sort -k6bnr > /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/10_richest_countries.txt
## extract the 10 poorest countries from the top of the file and present them in ascending order, i.e., the country with the lowest GDP on top
head -10 /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/2007_sortedbyGDP.txt > /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/10_poorest_countries.txt
## merge the two results in one file, starting with the richest countries and a file with two empty lines
cat /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/10_richest_countries.txt /mnt/c/users/sreev/KEYS/KEYS_Training/Scripts/EmptyLine.txt /mnt/c/users/sreev/KEYS/KEYS_Training/Scripts/EmptyLine.txt /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/10_poorest_countries.txt > /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/RichCountries_PoorCountries.txt
## add a header line with column headers
cat /mnt/c/users/sreev/KEYS/KEYS_Training/Scripts/EmptyLine.txt  /mnt/c/users/sreev/KEYS/KEYS_Training/Data/Data_Workfolder/gapminder_by_country/header.txt /mnt/c/users/sreev/KEYS/KEYS_Training/Scripts/EmptyLine.txt  /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/RichCountries_PoorCountries.txt > /mnt/c/users/sreev/KEYS/KEYS_Training/Analyses/RichCountries_PoorCountries_h.txt
## END SCRIPT     END SCRIPT     END SCRIPT     END SCRIPT     END SCRIPT     END SCRIPT     