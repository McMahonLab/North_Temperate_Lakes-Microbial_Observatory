::Script to run UniFrac in one go

::Some commands to avoid overwriting system variables
SETLOCAL ENABLEEXTENSIONS
SET me=%~n0
SET parent=%~dp0

::Run script to make count file and parse fasta file to remove sequences with no reads in these samples
"C:\Program Files\R\R-3.0.2\bin\x64\Rscript.exe" C:\Users\amlinz16\Desktop\North_Temperate_Lakes-Microbial_Observatory\UniFrac_analysis\02OTUtable_to_countfile.R %1 %2 %3

::Test to see if I can run mothur commands in the command prompt
C:\Users\amlinz16\Desktop\Programs\Mothur.win_64\mothur\mothur.exe C:\Users\amlinz16\Desktop\North_Temperate_Lakes-Microbial_Observatory\UniFrac_analysis\01unifrac_mothur_commands.txt


@EXIT /B 0