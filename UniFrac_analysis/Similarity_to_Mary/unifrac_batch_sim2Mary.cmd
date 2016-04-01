::Script to run UniFrac in one go

::Some commands to avoid overwriting system variables
SETLOCAL ENABLEEXTENSIONS
SET me=%~n0
SET parent=%~dp0

::Run script to make count file and parse fasta file to remove sequences with no reads in these samples
"C:\Program Files\R\R-3.2.4revised\bin\x64\Rscript.exe" C:\Users\amlinz16\Desktop\North_Temperate_Lakes-Microbial_Observatory\UniFrac_analysis\Similarity_to_Mary\01OTUtable_to_countfile_sim2Mary.R %1 %2 %3

::Run mothur commands
C:\Users\amlinz16\Desktop\Programs\Mothur.win_64\mothur\mothur.exe C:\Users\amlinz16\Desktop\North_Temperate_Lakes-Microbial_Observatory\UniFrac_analysis\Similarity_to_Mary\02unifrac_mothur_commands.txt

::Rename output
ren temp.trewsummary %1-%2-%3.trewsummary

::Get average UniFrac distance for each date in the mixed lake
"C:\Program Files\R\R-3.2.4revised\bin\x64\Rscript.exe" C:\Users\amlinz16\Desktop\North_Temperate_Lakes-Microbial_Observatory\UniFrac_analysis\03parse_unifrac_output_sim2Mary.R %1 %2 %3

@EXIT /B 0