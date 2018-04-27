try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os

packageName = "pyresid"
packageDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          packageName)

__version__ = "0.5"

setup(# package information
      name=packageName,
      version=__version__,
      author="Rob Firth",
      author_email="robert.firth@stfc.ac.uk",
      url="",
      description="Python tools for mining Protein Residuals from Fulltext articles using PMC number, ePMC and PDB",
      long_description='''Python tools for mining Protein Residuals from Fulltext articles using PMC number, ePMC and PDB.
  
                        BSD 3-Clause License
                        
                        Copyright (c) 2018, Robert Elliot Firth (Science and Techology Faclities Council)
                        All rights reserved.
                        
                        Redistribution and use in source and binary forms, with or without
                        modification, are permitted provided that the following conditions are met:
                        
                        * Redistributions of source code must retain the above copyright notice, this
                          list of conditions and the following disclaimer.
                        
                        * Redistributions in binary form must reproduce the above copyright notice,
                          this list of conditions and the following disclaimer in the documentation
                          and/or other materials provided with the distribution.
                        
                        * Neither the name of the copyright holder nor the names of its
                          contributors may be used to endorse or promote products derived from
                          this software without specific prior written permission.
                        
                        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
                        AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
                        IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
                        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
                        FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
                        DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
                        SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
                        CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
                        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
                        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.''',

      packages=[packageName,],
      package_dir={packageName: packageName},
      package_data={packageName:['mmCIF/*', "PDB/*"]},
      install_requires=["bs4", "scipy", "numpy", "matplotlib", "spacy", "lxml", "biopython", "pycifrw",
                        "fuzzywuzzy", "python-Levenshtein"]#,
    )
