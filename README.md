# CSHL

## Installation

### Mac

  1. Install [MacPorts](https://macports.org)
  2. Install required programs and libraries:
    
     ```bash
     sudo port install stimfit py27-stfio py27-scipy py27-matplotlib
     ```
  3. Get the script:
    
     ```bash
     cd
     git clone https://github.com/neurodroid/CSHL.git
     ```
  4. Place a symbolic link to the script in Python's library directory:
    
     ```bash
     sudo ln -s ~/CSHL/cshl.py /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/
     ```
     
  5. Whenever there's been a change to the script, update it with git:
  
     ```bash
     cd ~/CSHL
     git pull
     ```

### Windows

  1. Install the latest version of [stimfit](https://github.com/neurodroid/stimfit/releases).
  
  2. [Download the script](https://github.com/neurodroid/CSHL/raw/master/cshl.py) directly to `C:\Python27\Lib\site-packages`
   
  3. Repeat step 2 whenever there's been a change to the script
