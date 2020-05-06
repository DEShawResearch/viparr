Import('env')
import os

if True:
    cpp=[]
    flg=[]
    for p in env['CPPPATH']:
        if p.startswith('/proj') or p.startswith('/gdn'):
            flg.append('-I%s' % p)
        else:
            cpp.append(p)
    env.Replace(CPPPATH=cpp)
    env.Append(CFLAGS=flg, CXXFLAGS=flg)

env.Append(
    CPPDEFINES=['BOOST_FILESYSTEM_VERSION=3', 'BOOST_SYSTEM_NO_DEPRECATED'],
    CCFLAGS=['-O2', '-Wall', '-g', '-std=c++11', '-Wno-unused-local-typedefs'],
    LIBS=['boost_filesystem', 'boost_system', 'boost_program_options', 'msys', 'msys-core'],
    )

for d in 'src', 'python', 'tools':
    env.SConscript('%s/SConscript' % d)

env.AddShare('modules.txt')
env.AddShare('env.sh')

if 'BUILD_WHEEL' in os.environ:
    env['WHEEL_DIR'] = 'wheel'
    env.AddWheel('pyproject.toml', pyver='37')
