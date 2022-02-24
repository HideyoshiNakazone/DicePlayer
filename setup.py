import PyInstaller.__main__

name = 'diceplayer'

PyInstaller.__main__.run([
    'diceplayer/__main__.py',
    '--onefile',
    '-n{}'.format(name)
])