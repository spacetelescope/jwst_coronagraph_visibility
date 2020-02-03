# -*- mode: python -*-
import sys
from PyInstaller.utils import hooks

block_cipher = None

hiddenimports = ["tkinter", "tkinter.filedialog"]

a = Analysis(['run_jwst_coronagraph_visibility.py'],
             pathex=['.'],
             binaries=None,
             datas=hooks.collect_data_files('pysiaf'),
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
          name='run_jwst_coronagraph_visibility',
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
               name='run_jwst_coronagraph_visibility')

app = BUNDLE(coll,
             name='JWST Coronagraph Visibility Tool.app',
             icon='./jwst.icns',
             bundle_identifier='edu.stsci.jwst_coronagraph_visibility',
             info_plist={'NSHighResolutionCapable': 'True', 'CFBundleShortVersionString':'0.4.2'},
             )
