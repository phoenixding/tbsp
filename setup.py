from setuptools import setup

setup(  name='tbsp',
		version='1.0.0',
		description='SNP-based trajectory inference',
		author='Jun Ding',
		author_email='jund@andrew.cmu.edu',
		url="https://github.com/phoenixding/tbsp",
		license='MIT',
		packages=['tbsp'],
		entry_points={'console_scripts':['tbsp=tbsp.tbsp:main']},
		install_requires=['scipy','numpy','scikit-learn','matplotlib>=2.2.3','networkx>=2.2','Biopython>=1.72','pyBigWig'],
		classifiers=[
			'License :: OSI Approved :: MIT License',
			'Programming Language :: Python :: 2',
			'Programming Language :: Python :: 3',
		],
		)
		
