from setuptools import setup

setup(

	name="cprisma",
	version="1.0",
	description="This program gives color to multiple sequence alignment based on an input of numerical data.",
	author="Sergio Alejandro Poveda Cuevas",
	author_email="seapovedac@gmail.com",
	url="https://github.com/seapovedac/cprisma",
	license="GPLv3",
	packages=["cprisma"],
	include_package_data=True,
	package_data={'cprisma': ['*.csv']},
	entry_points={'console_scripts': ['cprisma = cprisma.main_run:main', ]},
	python_requires='>=3.7'

)
