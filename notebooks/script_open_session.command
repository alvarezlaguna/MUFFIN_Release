#! /bin/bash

cd /Users/alvarez/Documents/1_Simulation/DSMC_0D_ElectricField

# Open the folder
open ..
#open -a Terminal "../"
# Open VScode
open -a /Applications/Visual\ Studio\ Code.app "/Users/alvarez/Documents/1_Simulation/MUFFIN"

#open -a Terminal "`pwd`"

# Open one tab for the notebook
open -a Terminal "../"

osascript -e 'tell application "Terminal" to do script "cd /Users/alvarez/Documents/1_Simulation/MUFFIN/build; conda activate muffin-env;" in tab 1 of front window'

# Open one tab with Results
osascript -e 'tell application "Terminal" to activate' -e 'tell application "System Events" to tell process "Terminal" to keystroke "t" using command down'               
osascript -e 'tell application "Terminal" to do script "cd /Users/alvarez/Documents/1_Simulation/MUFFIN/Results; conda activate muffin-env; " in tab 1 of front window'

# Open one tab to Jupyter
osascript -e 'tell application "Terminal" to activate' -e 'tell application "System Events" to tell process "Terminal" to keystroke "t" using command down'         
osascript -e 'tell application "Terminal" to do script "conda activate muffin-env; jupyter-notebook" in tab 1 of front window'


# Open Firefox
open -n /Applications/Firefox.app 
open -a /Applications/Firefox.app https://www.overleaf.com/project/636cc53ea185d52d18768efe
open -a /Applications/Firefox.app https://github.com/LaboratoryOfPlasmaPhysics/MUFFIN

