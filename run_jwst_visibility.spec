# -*- mode: python -*-
import sys

block_cipher = None

hiddenimports = ['tkinter',]
if sys.version_info[0] == 2:
    hiddenimports.append('FileDialog')
else:
    hiddenimports.append('tkinter.filedialog')

a = Analysis(['run_jwst_visibility.py'],
             pathex=['.'],
             binaries=None,
             hiddenimports=hiddenimports,
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='run_jwst_visibility',
          debug=False,
          strip=False,
          upx=True,
          console=False )

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='run_jwst_visibility')

app = BUNDLE(coll,
             name='JWST Visibility Calculator.app',
             icon='./jwst.icns',
             bundle_identifier='edu.stsci.jwst_visibility',
             info_plist={'NSHighResolutionCapable': 'True'},
             )
