# Getting to grips with working in terminal 

A lot of bioinformatics pipelines and softwares are run directly in terminal, familiarity with bash scripting will make learning these softwares hopefully a little easier 

A great resource, as always, is the Pyhton For Biologist webpage from Helfrid Hochegger.

Here is a link to the page specifically about bash:
https://python-for-biologists.vercel.app/dev-tools/cursor-setup/terminal

Another excellent bash basics resource is here: 
https://media.datacamp.com/legacy/image/upload/v1700047731/Marketing/Blog/Bash_Cheat_Sheet.pdf

---
## Setting up your computer (Mac or Linux)

These instructions are also copied and adapted from the Hochegger PFB webpage: https://python-for-biologists.vercel.app/dev-tools/cursor-setup#:~:text=Complete%20Setup%20Instructions

### Overview of install
1. ✅ Warp Terminal (modern terminal with AI features)
2. ✅ Git (for downloading materials)
   Please only proceed to the next step if you do not already have python and/or a package manager installed
3. ✅ mambaforge install to manage packages

---
#### Step 1: Warp Terminal install

Warp is a modern terminal with AI features that makes command-line work much easier for beginners. It also allows you to more easily navigate without using the arrow pad compared to the normal terminal.

**For Macs:**
1. Visit https://www.warp.dev/
2. Click “Download for Mac”
3. Open the downloaded .dmg file
4. Drag Warp to your Applications folder 5. Open Warp from Applications
6. Sign in with your email (free account)

---
#### Step 2: Checking for Git
Before installing anything, let’s check if you already have Git.

Open Warp (or terminal) and run:
```bash
git --version
```

If you see a version number (like git version 2.39.0):
• ✅ You’re good to go! Skip to Step 3. If you see “command not found”:
• Continue with the installation instructions below.

First, install **Homebrew**
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/
        Homebrew/install/HEAD/install.sh)"
```
Follow the on-screen instructions and enter your password when prompted. 

Then install **Git**:
```bash
brew install git
```

Check your install

```bash
git --version
```
---

#### Step 3:
Installing conda/mamba forge

Follow the instructions from this webpage: https://github.com/conda-forge/miniforge 
to install miniforge

After this is done, we should be ready to set up conda/mamba environments to work in. This is overall more stable and better long term compared to using something like Anaconda (also a package manager)

**Critical point**
Please remember that anytime you decide from now on to install a package, please do so within a python environment, and use `conda install ` whenever you can.
If you plan to use an IDE like Spyder, do NOT install from the Spyder console, install instead from the terminal/Warp.

This should help to keep everything nice and tidy, and hopefully prevent your environment from crashing :)

---

## Bash worksheet:
To help you get used to some of the commands you may need with bash scripting, please go through the **L3_bash_worksheet.ipynb** in this repository

---


