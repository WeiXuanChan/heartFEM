from setuptools import setup
setup(
  name = 'heartFEM',         # How you named your package folder (MyLib)
  packages = ['heartFEM','heartFEM.lcleeHeart'],  
  include_package_data=True,   
  package_data={
        'ngspice':['*'],
        'lcleeHeart':['*','vtk_py/*']
        },
  version = '3.5.1',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Optimization of patient specific heart model with myocardium FEM and cardiovascular network Windkessel.',   # Give a short description about your library
  author = 'Wei Xuan Chan',                   # Type in your name
  author_email = 'w.x.chan1986@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/WeiXuanChan/heartFEM',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/WeiXuanChan/heartFEM/archive/v3.5.1.tar.gz',    
  keywords = ['medical', 'cardiac'],   # Keywords that define your package best
  install_requires=['numpy','matplotlib','scipy'],
  extras_require = {
        'fem analysis software':  ['fenics'],
        'circuit analysis':  ['ngspice'],
        'data storage and processing':  ['vtk']
  },
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package    
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',    
    'License :: OSI Approved :: MIT License',   # Again, pick a license    
    'Programming Language :: Python :: 3.6',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.7',
  ],
)
