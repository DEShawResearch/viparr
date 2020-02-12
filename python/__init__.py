from __future__ import print_function

import msys
from . import _viparr
from .version import version as __version__
import os
import atexit
import warnings

def _comb_rule_to_func(key, _VDWCombRule):
    def f(vdw_params_A, vdw_params_B, lj_scale):
        return _VDWCombRule(vdw_params_A, vdw_params_B, lj_scale)
    f.__name__ = key
    return f

def get_ffpath():
    """Return the forcefield search path"""
    return os.getenv("VIPARR_FFPATH", "")

def list_forcefields(ffpath):
    """Return a dict mapping forcefield name to normalized full path

    Notes
    -----
        A warning is issued if two forcefields with the same name are found in multiple locations
    """
    ffs = dict()
    seen = set()
    for ffdir in ffpath.split(':'):
        realdir = os.path.realpath(ffdir)
        if realdir in seen: continue
        seen.add(realdir)
        for name in os.listdir(realdir):
            if name in ffs:
                warnings.warn("Forcefield '%s' appears in two different locations in your VIPARR_FFPATH: '%s' and '%s'" % (
                    name, realdir, ffs[name]))
            rules_file = '%s/%s/rules' % (realdir, name)
            if os.path.exists(rules_file):
                ffs[name] = '%s/%s' % (realdir, name)
    return ffs

def find_forcefield(name):
    """Return the full path to the forcefield with the given name, or raise a RuntimeError.
    """
    if os.path.isdir(name):
        return os.path.abspath(name)
    ffpath = get_ffpath()
    if not ffpath:
        raise RuntimeError("Could not find forcefield '%s' because VIPARR_FFPATH is unset; did you load a viparr-ff/data module?" % name)
    ffdir = list_forcefields(ffpath).get(name)
    if not ffdir:
        raise RuntimeError("Forcefield '%s' not found; search path was '%s'" % (name, ffpath))
    return ffdir


########################## Forcefield classes ##################################
class Forcefield(object):
    """Representation of a forcefield directory.

    A :class:`Forcefield` object consists of:

    (1) a :class:`Rules` object, containing the nonbonded info, exclusion and
        scaled-pairs info, and list of forcefield plugins to apply
    
    (2) a :class:`TemplateTyper`, :class:`SmartsTyper`, or
        :class:`ScoredSmartsTyper` object, containing information to assign atom
        types, extra exclusions, impropers, cmaps, and/or pseudo particles to
        the system
    
    (3) a map from param table name to list of :class:`msys.Param` objects,
        which are pointers to rows of a collection of static shared parameter
        tables, indicating which shared parameter tables and which parameters in
        those tables belong to this forcefield
    
    (4) an optional list of (non-shared) cmap parameter tables
 
    The :class:`Forcefield` class contains a static dictionary of shared
    parameter tables, with static accessors and modifiers. :class:`viparr` is
    designed such that all parameters from forcefields and systems should
    belong to this static dictionary, and all parameters of the same type (e.g.
    stretch-harm) should belong to a single table in this dictionary (e.g. with 
    key 'stretch-harm'). When :func:`ImportForcefield` or :func:`ImportParams` 
    is called, new parameters from the imported forcefield are automatically 
    merged into this static dictionary. If a chemical system containing existing
    forcefield information is loaded using :class:`msys`, its parameters may be
    merged into this static dictionary using the global :func:`AddSystemTables`
    function.

    The :class:`Forcefield` class also maintains a static registry of supported
    parameter-matching plugins.
    """

    @classmethod
    def _from_boost(cls, _Forcefield):
        if _Forcefield is None:
            return None
        ff = cls(Rules(), TemplateTyper())
        ff._Forcefield = _Forcefield
        return ff

    @staticmethod
    def HasParamTable(name):
        """Whether a given param table is in the static param table dictionary.

        Arguments:
            name -- str

        Returns: bool
        """
        return _viparr.Forcefield.HasParamTable(name)

    @staticmethod
    def ParamTable(name):
        """Get a static param table by name.

        Arguments:
            name -- str

        Returns: :class:`msys.ParamTable`
        """
        return msys.ParamTable(_viparr.Forcefield.ParamTable(name))

    @staticmethod
    def AddParamTable(name, table):
        """Add a param table to the static param table dictionary.

        Arguments:
            name -- str

            table -- :class:`msys.ParamTable`

        """
        _viparr.Forcefield.AddParamTable(name, table._ptr)

    @staticmethod
    def AllParamTables():
        """List of keys (table names) of the static param table dictionary.

        Returns: [str, ..., str]
        """
        return _viparr.Forcefield.AllParamTables()

    @staticmethod
    def ClearParamTables():
        """Clear the static param table dictionary."""
        _viparr.Forcefield.ClearParamTables()

    class Plugin:
        """A wrapper class around forcefield plugin functions.

        A typical forcefield plugin function matches a particular type of term
        (e.g. bonds or angles) in a :class:`TemplatedSystem` to parameters in a
        particular param table of the forcefield, creating an
        :class:`msys.TermTable` in the system to save the parametrized terms.
        At a minimum, the :class:`Plugin` object wraps such a function

            apply(:class:`TemplatedSystem`, :class:`Forcefield`) -> None.

        Adding a :class:`Plugin` object that wraps such an apply function to
        the static plugin registry of the :class:`Forcefield` class is all that
        is required for use of this new plugin in viparr execution.

        A :class:`Plugin` can optionally specify a list of prerequisite plugins
        which, if also present in the :class:`Rules` of a forcefield, must be
        applied before this one (e.g. 'angles' must be applied before
        'inplanewags'). :class:`viparr` will check that these prerequirement
        conditions are satisfied during its execution.
        """
        @classmethod
        def _from_boost(cls, _Plugin):
            if _Plugin is None:
                return None
            plugin = cls(None,None)
            plugin._Plugin = _Plugin
            return plugin

        def __init__(self, matchf, compilef=lambda x:None, prerequisites=[]):
            """Create from apply function.

            Arguments:
                apply -- f(:class:`TemplatedSystem`, :class:`Forcefield`) -> None

                prerequisites -- [str, ..., str], list of plugin names

            """
            def apply_c(tsystem, ff):
                matchf(TemplatedSystem._from_boost(tsystem),
                       Forcefield._from_boost(ff))
            def compile_c(sys):
                compilef(sys)
            self._Plugin = _viparr.Forcefield.Plugin(apply_c,compile_c)
            self._Plugin.setPrerequisites(prerequisites)

        @property
        def prerequisites(self):
            """ List of plugin prerequisites (plugin names)"""
            return self._Plugin.getPrerequisites()

        @prerequisites.setter
        def prerequisites(self, prerequisites):
            self._Plugin.setPrerequisites(prerequisites)

        def apply(self, tsystem, ff):
            """Apply the contained plugin function

            Arguments:
                tsystem -- :class:`TemplatedSystem`

                ff -- :class:`Forcefield`

            """
            self._Plugin.match(tsystem._TemplatedSystem, ff._Forcefield)

        def compile(self, system):
            self._Plugin.compile(system._ptr)

    @classmethod
    def GetPluginRegistry(cls):
        """Return the registry of all known plugins, indexed by name.

        When the :class:`viparr` library is imported, the registry is populated
        by default with a list of viparr's built-in plugins, including 'angles',
        'bonds', 'cmap', 'exclusions', 'impropers', 'mass', 'propers',
        'pseudopol_fermi', 'ureybradley', 'vdw1', and 'virtuals'.

        Returns: { str: :class:`Plugin`, ..., str: :class:`Plugin` }
        """
        return dict((key, cls.Plugin._from_boost(_Plugin)) for (key,
            _Plugin) in _viparr.Forcefield.PluginRegistry.items())

    @staticmethod
    def AddPlugin(key, plugin):
        """Add a :class:`Plugin` to the registry.

        Arguments:
            key -- str

            plugin -- :class:`Plugin`

        """
        registry = _viparr.Forcefield.PluginRegistry
        registry[key] = plugin._Plugin
        _viparr.Forcefield.PluginRegistry = registry

    @staticmethod
    def DelPlugin(key):
        """Delete a :class:`Plugin` from the registry.

        Arguments:
            key -- str

        """
        registry = _viparr.Forcefield.PluginRegistry
        del registry[key]
        _viparr.Forcefield.PluginRegistry = registry

    def __init__(self, rules, typer):
        """Construct with no param tables and no cmap tables.

        Arguments:
            rules -- :class:`Rules`

            typer -- :class:`TemplateTyper` or :class:`SmartsTyper` or :class:`ScoredSmartsTyper`

        """
        self._Forcefield = _viparr.Forcefield(rules._Rules, typer._Typer)

    def copy(self):
        """Create a copy of this forcefield.

        The forcefield copy shares this forcefield's :class:`Rules` and
        :class:`Typer` objects but has its own lists of :class:`msys.Param`
        objects (pointers into rows of the shared param tables).

        Returns: :class:`Forcefield`
        """
        _ff = _viparr.Forcefield.Copy(self._Forcefield)
        return Forcefield._from_boost(_ff)

    def __eq__(self, other):
        try:
            return self._Forcefield == other._Forcefield
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._Forcefield.__hash__()

    @property
    def name(self):
        """The name of the forcefield.
        
        Defaults to the system path for forcefields created by
        :func:`ImportForcefield`.
        """
        return self._Forcefield.name

    @name.setter
    def name(self, val):
        self._Forcefield.name = val

    @property
    def rules(self):
        """The contained :class:`Rules` object."""
        return Rules._from_boost(self._Forcefield.rules())

    @rules.setter
    def rules(self, val):
        self._Forcefield.resetRules(val)

    @property
    def typer(self):
        """The contained :class:`TemplateTyper` object.
        
        May be a :class:`SmartsTyper` or a :class:`ScoredSmartsTyper`."""
        return TemplateTyper._from_boost(self._Forcefield.typer())
    
    @typer.setter
    def typer(self, val):
        self._Forcefield.resetTyper(val._Typer)

    @property
    def paramTables(self):
        """List of names of param tables belonging to this forcefield."""
        return self._Forcefield.paramTables()

    def params(self, name):
        """Retrieve params belonging to this forcefield for a given table name.

        If 'name' does not correspond to a param table in this forcefield,
        return value is [].

        Arguments:
            name -- str

        Returns: [:class:`msys.Param`, ..., :class:`msys.Param`]
        """
        rowIDs = self._Forcefield.rowIDs(name)
        table = Forcefield.ParamTable(name)
        return [msys.Param(table._ptr, id) for id in rowIDs]

    def findParams(self, name, **kwds):
        """Find parameters from a given table with the given property values.

        Names of keyword arguments should correspond to properties of the
        param table with the given name; this function will return all params
        matching all of the given values in this forcefield from this table.

        Returns: [:class:`msys.Param`, ..., :class:`msys.Param`]
        """
        if not kwds:
            return []
        rowIDs = self._Forcefield.rowIDs(name)
        if len(rowIDs) == 0:
            return []
        for k, v in kwds.items():
            if type(v) == int:
                rowIDs = _viparr.Forcefield.FilterParamsInt(name, rowIDs, k, v)
            elif type(v) == float:
                rowIDs = _viparr.Forcefield.FilterParamsFloat(name, rowIDs, k,
                        v)
            elif type(v) == str:
                rowIDs = _viparr.Forcefield.FilterParamsString(name, rowIDs, k,
                        v)
            else:
                raise RuntimeError("Type of value for keyword '" + k + \
                        "' must be int, float, or str")
            if len(rowIDs) == 0:
                return []
        table = Forcefield.ParamTable(name)
        return [msys.Param(table._ptr, id) for id in rowIDs]

    def appendParam(self, name, param=None, **kwds):
        """Append a new parameter to a given table of this forcefield.

        If 'param' is given, it should belong to the param table of the
        given name, and this existing param is added. Otherwise, a new parameter
        will be added to the global table with values specified by the given
        keyword arguments, and this parameter will be added to the forcefield.

        Arguments:
            name -- str

            param -- :class:`msys.Param`

        Returns: msys.Param
        """
        if param is None:
            if not Forcefield.HasParamTable(name):
                Forcefield.AddParamTable(name, msys.CreateParamTable())
                for k, v in kwds.items():
                    Forcefield.ParamTable(name).addProp(k, type(v))
            param = Forcefield.ParamTable(name).addParam(**kwds)
        if not Forcefield.HasParamTable(name) or \
                param._ptr != self.ParamTable(name)._ptr:
            raise RuntimeError("Parameter " + repr(param) + \
                    " does not belong to this table")
        self._Forcefield.appendParams(name, [param.id])
        return param

    def delParams(self, name, params=None, **kwds):
        """Delete parameters from a given table for this forcefield.

        If 'params' are given, they will be removed from this forcefield.
        Otherwise, findParams will be used to find all parameters with the
        properties and values given by the keyword arguments, and all such
        parameters will be removed.

        Arguments:
            name -- str

            param -- :class:`msys.Param` or [:class:`msys.Param`, ... ]

        """
        if params is None:
            params = self.findParams(name, **kwds)
        elif type(params) == msys.Param:
            params = [params]
        for p in params:
            if p._ptr != Forcefield.ParamTable(name)._ptr:
                raise RuntimeError("Parameter " + p.__repr__() + \
                        " does not belong to this table")
        self._Forcefield.delParams(name, [p.id for p in params])

    def replaceParam(self, name, old_param, new_param=None, **kwds):
        """Replace a parameter with another parameter in the same position.

        'old_param' should be a parameter in table 'name' of this forcefield.
        If 'new_param' is given, it should also be a parameter in table 'name'
        of this forcefield, and it will replace 'old_param'. Otherwise, a new
        parameter will be appended to the global table using values given by
        the keyword arguments, and this will replace 'old_param'.

        Arguments:
            name -- str

            old_param -- :class:`msys.Param`

            new_param -- :class:`msys.Param`

        Returns: :class:`msys.Param`
        """
        if old_param._ptr != Forcefield.ParamTable(name)._ptr:
            raise RuntimeError("Parameter " + old_param.__repr__() + \
                    " does not belong to this table")
        if new_param is None:
            new_param = Forcefield.ParamTable(name).addParam(**kwds)
        if new_param._ptr != Forcefield.ParamTable(name)._ptr:
            raise RuntimeError("Parameter " + new_param.__repr__() + \
                    " does not belong to this table")
        self._Forcefield.replaceParam(name, old_param.id, new_param.id)
        return new_param

    @property
    def cmap_tables(self):
        """All contained cmap parameter tables."""
        tables = self._Forcefield.cmapTables()
        return [msys.ParamTable(table) for table in tables]

    def addCmapTable(self, table):
        """Add a new cmap parameter table."""
        self._Forcefield.addCmapTable(table._ptr)

    def __str__(self):
        output = 'Name: ' + self.name + '\n\n'
        output += 'Rules:\n' + self.rules.__str__()
        output += 'Typer:\n' + self.typer.__str__()
        output += 'Param tables:\n'
        for table in self.paramTables:
            output += table + '\n'
            output += PrintParams(self.params(table)) + '\n'
        if len(self.cmap_tables) > 0:
            output += 'Cmap tables:\n'
            for ind, table in enumerate(self.cmap_tables):
                output += 'Table ' + str(ind) + ':\n'
                output += PrintParams(table.params) + '\n'
        return output

    def __repr__(self):
        return "<Forcefield '%s'>" % self.name

# Viparr's static Forcefield registry can contain Python-created objects
# which need to be fully destructed before Python shuts down.  Install
# a handler to make sure that happens.
atexit.register(Forcefield.ClearParamTables)
atexit.register(_viparr.Forcefield.ClearPlugins)

# Do the same for the vdw func and rules registry.
atexit.register(_viparr.Rules.ClearVdwFuncRegistry)
atexit.register(_viparr.Rules.ClearVdwCombRuleRegistry)

class TemplatedSystem(object):
    """A wrapper around :class:`msys.System` to hold template information.
    
    Serves as a vehicle for communication between the atom typing and parameter
    matching steps of :class:`viparr` execution. Can store atom types and bond
    aromaticity, as well as lists of typed atoms and bonds, angles, dihedrals,
    exclusions, impropers, cmaps, and pseudos involving those typed atoms that
    are needed for parameter matching. :class:`viparr` execution operates by
    applying a forcefield to a system to generate a :class:`TemplatedSystem`
    containing information for a subset of typed fragments of the system, using
    this :class:`TemplatedSystem` to match the parameters of the forcefield and
    apply forcefield plugins on these typed fragments, and then discarding this
    :class:`TemplatedSystem` before applying the next forcefield to parametrize
    the remaining fragments.
    
    :class:`TemplatedSystem` objects are also used to represent forcefield
    templates, storing atom-type information, impropers, extra exclusions, etc. 
    as defined in the templates.
    
    Pseudo particles (including virtuals and drudes) have atomic number 0.
    When used to represent Forcefield templates, externally bonded atoms have
    atomic number -1.
    """

    @classmethod
    def _from_boost(cls, _TemplatedSystem):
        if _TemplatedSystem is None:
            return None
        tsystem = cls()
        tsystem._TemplatedSystem = _TemplatedSystem
        return tsystem

    def __init__(self, system=None):
        """ Construct around an :class:`msys.System` object.

        The contained :class:`msys.System` object is stored by reference; one
        may construct multiple :class:`TemplatedSystem` wrappers around the same
        :class:`msys.System`. If 'system' is None, constructs the wrapper around
        a new, empty system.

        Arguments:
            system -- :class:`msys.System`

        """
        if system:
            self._TemplatedSystem = _viparr.TemplatedSystem(system._ptr)
        else:
            self._TemplatedSystem = _viparr.TemplatedSystem()
    
    def __eq__(self, other):
        try:
            return self._TemplatedSystem == other._TemplatedSystem
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._TemplatedSystem.__hash__()

    @property
    def system(self):
        """The contained :class:`msys.System` object."""
        return msys.System(self._TemplatedSystem.system())

    def clone(self, seltext=None):
        """Copy the viparr.TemplatedSystem, returning the copy.

        If seltext is provided, it should be a valid VMD atom
        selection, and only the selected atoms will be cloned. If
        seltext is a sequence (and not a string), it will be treated
        as a list of atom ids.
        
        Arguments:
            seltext -- A sequence of atom ids or a VMD atom selection

        Returns:
            :class:`viparr.TemplatedSystem`
            
        """
        if seltext is None:
            ids = [atom.id for atom in self.system.atoms]
        elif isinstance(seltext, str):
            ids = self.system.selectIds(seltext)
        else:
            ids = []         # formerly _msys.IdList()
            for a in seltext: ids.append(a)

        tsystem = self.__class__()
        tsystem._TemplatedSystem = self._TemplatedSystem.clone(ids)
        return tsystem

    @property
    def name(self):
        """Name of template, equal to name of first residue in the system."""
        return self.system.residue(0).name

    @name.setter
    def name(self, val):
        self.system.residue(0).name = val

    def btype(self, atom):
        """Return the bonded atom type of an atom.

        Arguments:
            atom -- :class:`msys.Atom`

        Returns: str

        """
        return self._TemplatedSystem.btype(atom._id)

    def nbtype(self, atom):
        """Return the nonbonded atom type of an atom.

        This is generally the same as btype for SMARTS based atomtyping but may
        be different in forcefield templates.

        Arguments:
            atom -- :class:`msys.Atom`

        Returns: str

        """
        return self._TemplatedSystem.nbtype(atom._id)

    def pset(self, atom):
        """Return the pset of an atom.

        Return value is generally '' if the atom is not a pseudo particle.

        Arguments:
            atom -- :class:`msys.Atom`

        Returns: str

        """
        return self._TemplatedSystem.pset(atom._id)

    def setTypes(self, atom, btype, nbtype, pset=''):
        """Set the bonded type, nonbonded type, and pset of an atom.

        'pset' should generally be set to '' if the atom is not a pseudo
        particle.

        Arguments:
            atom -- :class:`msys.Atom`

            btype -- str

            nbtype -- str

            pset -- str

        """
        self._TemplatedSystem.setTypes(atom._id, btype, nbtype, pset)

    @property
    def hash(self):
        """:class:`msys.Graph` hash of all non-pseudo atoms in the system."""
        return self._TemplatedSystem.hash()

    @property
    def graph(self):
        """:class:`msys.Graph` of all non-pseudo atoms in the system."""
        return msys.Graph._from_boost(self._TemplatedSystem.graph(),
                self._TemplatedSystem.system())

    def addTypedAtom(self, atom):
        """Add an atom to the list of typed atoms.

        This does not add any atoms to the contained :class:`msys.System`; it
        only adds an atom already existing in the :class:`msys.System` to the
        :class:`TemplatedSystem`'s list of typed atoms (for which charge, vdw,
        etc. parameters will be matched).

        Arguments:
            atom -- :class:`msys.Atom`

        """
        self._TemplatedSystem.addTypedAtom(atom._id)

    def addNonPseudoBond(self, atoms):
        """Add a pair of atoms to the list of non-pseudo bonds.

        This does not add any bonds to the contained :class:`msys.System`; it
        only adds a bond already existing in the :class:`msys.System` to the
        :class:`TemplatedSystem`'s list of bonds (for which stretch_harm, etc.
        parameters will be matched).

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.addNonPseudoBond([atom._id for atom in atoms])

    def addPseudoBond(self, atoms):
        """Add a pair of atoms to the list of pseudo bonds.

        This does not add any bonds to the contained :class:`msys.System`; it
        only adds a bond already existing in the :class:`msys.System` to the
        :class:`TemplatedSystem`'s list of pseudo bonds (used when matching
        pseudo parameters).  

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.addPseudoBond([atom._id for atom in atoms])

    def addAngle(self, atoms):
        """Add a triple of atoms to the list of angles.

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.addAngle([atom._id for atom in atoms])

    def addDihedral(self, atoms):
        """Add a quadruple of atoms to the list of dihedral angles.

        Arguments:
            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        """
        self._TemplatedSystem.addDihedral([atom._id for atom in atoms])

    def addExclusion(self, atoms):
        """Add a pair of atoms to the list of extra exclusions.

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.addExclusion([atom._id for atom in atoms])

    def addImproper(self, atoms):
        """Add a quadruple of atoms to the list of impropers.

        Arguments:
            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        """
        self._TemplatedSystem.addImproper([atom._id for atom in atoms])

    def addCmap(self, atoms):
        """Add an octuple of atoms to the list of cmaps.

        Arguments:
            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        """
        self._TemplatedSystem.addCmap([atom._id for atom in atoms])

    def addPseudoType(self, type, nsites):
        """Add the definition for a new type of pseudo particle.

        Arguments:
            type -- str

            nsites -- int, number of virtual sites including the pseudo itself

        """
        self._TemplatedSystem.addPseudoType(type, nsites)

    def addPseudoSites(self, type, atoms):
        """Add a pseudo atom and its virtual sites for a given pseudo type.

        This does not add any atoms to the contained msys.System; it only adds
        existing atoms to the TemplatedSystem's list of pseudo atom-tuples for
        the given pseudo type (which are used in matching pseudo parameters).

        Arguments:
            type -- str

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`], where
                atoms[0] is the pseudo itself and atoms[1],...atoms[n] are the
                site atoms

        """
        self._TemplatedSystem.addPseudoSites(type, [atom._id for atom in atoms])

    def removeTypedAtom(self, atom):
        """Remove an atom from the list of typed atoms.

        This does not remove any atoms to the contained :class:`msys.System`; it
        only removes an atom already existing in the :class:`msys.System` to the
        :class:`TemplatedSystem`'s list of typed atoms (for which charge, vdw,
        etc. parameters will be matched).

        Arguments:
            atom -- :class:`msys.Atom`

        """
        self._TemplatedSystem.removeTypedAtom(atom._id)

    def removeNonPseudoBond(self, atoms):
        """Remove a pair of atoms from the list of non-pseudo bonds.

        This does not remove any bonds to the contained :class:`msys.System`; it
        only removes a bond already existing in the :class:`msys.System` to the
        :class:`TemplatedSystem`'s list of bonds (for which stretch_harm, etc.
        parameters will be matched).

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.removeNonPseudoBond([atom._id for atom in atoms])

    def removePseudoBond(self, atoms):
        """Remove a pair of atoms from the list of pseudo bonds.

        This does not remove any bonds to the contained :class:`msys.System`; it
        only removes a bond already existing in the :class:`msys.System` to the
        :class:`TemplatedSystem`'s list of pseudo bonds (used when matching
        pseudo parameters).  

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.removePseudoBond([atom._id for atom in atoms])

    def removeAngle(self, atoms):
        """Remove a triple of atoms from the list of angles.

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.removeAngle([atom._id for atom in atoms])

    def removeDihedral(self, atoms):
        """Remove a quadruple of atoms from the list of dihedral angles.

        Arguments:
            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        """
        self._TemplatedSystem.removeDihedral([atom._id for atom in atoms])

    def removeExclusion(self, atoms):
        """Remove a pair of atoms from the list of extra exclusions.
        This method is insensitive to atom ordering:
        removeExclusion(a1, a2) has the same effect as removeExclusion(a2, a1).

        Arguments:
            atoms -- [:class:`msys.Atom`, :class:`msys.Atom`]

        """
        self._TemplatedSystem.removeExclusion([atom._id for atom in atoms])

    def removeImproper(self, atoms):
        """Remove a quadruple of atoms from the list of impropers.

        Arguments:
            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        """
        self._TemplatedSystem.removeImproper([atom._id for atom in atoms])

    def removeCmap(self, atoms):
        """Remove an octuple of atoms from the list of cmaps.

        Arguments:
            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        """
        self._TemplatedSystem.removeCmap([atom._id for atom in atoms])
        
    @property
    def typedAtoms(self):
        """List of typed atoms, including typed pseudo particles.

        Returns: [[:class:`msys.Atom`], ..., [:class:`msys.Atom`]]
        """
        idslist = self._TemplatedSystem.typedAtoms()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def nonPseudoBonds(self):
        """List of bonds between non-pseudo typed atoms.

        Returns: [[:class:`msys.Atom`, :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.nonPseudoBonds()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def pseudoBonds(self):
        """List of bonds between a pseudo and a non-pseudo typed atom.

        Returns: [[:class:`msys.Atom`, :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.pseudoBonds()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def angles(self):
        """List of angles between non-pseudo typed atoms.

        Returns: [[:class:`msys.Atom`, :class:`msys.Atom`, :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.angles()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def dihedrals(self):
        """List of dihedral angles between non-pseudo typed atoms.

        Returns: [[:class:`msys.Atom`, ..., :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.dihedrals()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def exclusions(self):
        """List of extra exclusions between pairs of typed atoms.

        Returns: [[:class:`msys.Atom`, :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.exclusions()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def impropers(self):
        """List of improper dihedral angles between non-pseudo typed atoms.

        Returns: [[:class:`msys.Atom`, ..., :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.impropers()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def cmaps(self):
        """List of cmap octuples for non-pseudo typed atoms.

        Returns: [[:class:`msys.Atom`, ..., :class:`msys.Atom`], ... ]
        """
        idslist = self._TemplatedSystem.cmaps()
        return [[msys.Atom(self.system._ptr, id) for id in ids] for ids in idslist]

    @property
    def pseudoTypes(self):
        """List of pseudo-type definitions and typed pseudo particles.

        Returns: [:class:`PseudoType`, ..., :class:`PseudoType`]
        """
        _pseudo_types = self._TemplatedSystem.pseudoTypes()
        return [TemplatedSystem.PseudoType._from_boost(_pseudo_type,
            self._TemplatedSystem) for _pseudo_type in _pseudo_types]

    class PseudoType(object):
        """Definition for a pseudo type and container for list of typed pseudos.
        """
        @classmethod
        def _from_boost(cls, _PseudoType, _TemplatedSystem):
            if _PseudoType is None:
                return None
            pseudo_type = cls()
            pseudo_type._PseudoType = _PseudoType
            pseudo_type._system = _TemplatedSystem.system()
            return pseudo_type

        def __repr__(self):
            return '<PseudoType %s>' % self._PseudoType.name

        @property
        def name(self):
            """Name of pseudo type."""
            return self._PseudoType.name

        @property
        def nsites(self):
            """Number of sites for pseudo type, including the pseudo itself."""
            return self._PseudoType.nsites

        @property
        def sites_list(self):
            """List of typed pseudos and their site atoms for this pseudo type.

            Returns: [[:class:`msys.Atom`, ... ], ...] where the first atom
                of each returned atom-tuple is the pseudo particle itself
                and the remaining atoms are the site atoms.

            """
            idslist = self._PseudoType.sites_list
            return [[msys.Atom(self._system, id) for id in ids] \
                    for ids in idslist]

class Rules(object):
    """Representation of a rules forcefield file.

    Contains the forcefield's nonbonded info, exclusion rule, ES and LJ
    scaling-factors for scaled pair interactions, list of plugins, and
    several boolean flags controlling viparr behavior.
    
    The :class:`Rules` class also contains static registries of supported VDW
    functions and VDW combine rules.
    """
    @classmethod
    def _from_boost(cls, _Rules):
        if _Rules is None:
            return None
        rules = cls()
        rules._Rules = _Rules
        return rules

    class VDWFunc(object):
        """Representation of a VDW functional form."""

        @classmethod
        def _from_boost(cls, _VDWFunc):
            if _VDWFunc is None:
                return None
            func = cls()
            func._VDWFunc = _VDWFunc
            return func

        def __init__(self):
            self._VDWFunc = _viparr.Rules.VDWFunc()

        def __eq__(self, other):
            try:
                return self._VDWFunc == other._VDWFunc
            except AttributeError:
                return False

        def __ne__(self, other):
            return not self.__eq__(other)

        def __hash__(self):
            return self._VDWFunc.__hash__()

        @property
        def vdw_table_name(self):
            """Name of the nonbonded table to be added to the system."""
            return self._VDWFunc.vdw_table_name

        @vdw_table_name.setter
        def vdw_table_name(self, val):
            self._VDWFunc.vdw_table_name = val

        @property
        def param_names(self):
            """Names of the parameters in the VDW param table."""
            return self._VDWFunc.param_names

        @param_names.setter
        def param_names(self, val):
            self._VDWFunc.param_names = val

        @property
        def pair_table_name(self):
            """Name of the nonbonded pairs table to be added to the system."""
            return self._VDWFunc.pair_table_name

        @pair_table_name.setter
        def pair_table_name(self, val):
            self._VDWFunc.pair_table_name = val

        @property
        def pair_param_names(self):
            """Names of the combined VDW pair parameters."""
            return self._VDWFunc.pair_param_names

        @pair_param_names.setter
        def pair_param_names(self, val):
            self._VDWFunc.pair_param_names = val

        @property
        def supported_rules(self):
            """List of all VDW combine rules supported by this functional form.

            Each combine rule should be a key in
            :func:`Rules.GetVDWCombRuleRegistry()`.
            """
            return self._VDWFunc.supported_rules

        @supported_rules.setter
        def supported_rules(self, val):
            self._VDWFunc.supported_rules = val


    @staticmethod
    def GetVDWFuncRegistry():
        """Return registry of supported VDW functional forms, indexed by name.

        The registry is populated by default with 'lj12_6_sig_epsilon' and
        'exp_6x'.

        Returns: {str: :class:`VDWFunc`, ..., str: :class:`VDWFunc`}
        """
        return dict([(key, Rules.VDWFunc._from_boost(func)) \
                for key, func in _viparr.Rules.VDWFuncRegistry.items()])

    @staticmethod
    def AddVDWFunc(key, vdw_func):
        """Add a new VDW functional form to the registry.

        This allows viparr to use a forcefield with the new functional form
        specified in its rules file.

        Arguments:
            key -- str

            vdw_func -- :class:`VDWFunc`

        """
        registry = _viparr.Rules.VDWFuncRegistry
        registry[key] = vdw_func._VDWFunc
        _viparr.Rules.VDWFuncRegistry = registry

    @staticmethod
    def DelVDWFunc(key):
        """Delete a VDW functional form from the registry.

        Arguments:
            key -- str

        """
        registry = _viparr.Rules.VDWFuncRegistry
        del registry[key]
        _viparr.Rules.VDWFuncRegistry = registry

    # We do not expose the VDWCombRule wrapper class in Python, but rather
    # convert default VDWCombRules to Python functions and provide an interface
    # to operate directly on Python functions
    _comb_rule_map = dict([(key, _comb_rule_to_func(key, rule)) \
            for key, rule in _viparr.Rules.VDWCombRuleRegistry.items()])

    @classmethod
    def GetVDWCombRuleRegistry(cls):
        """Return registry of supported VDW combine rules, indexed by name.

        A VDW combine rule is a function

            f(params_A, params_B, scale_factor) -> combined_params,

        where params_A and params_B are lists of VDW parameter values
        corresponding (in sequential order) to the properties
        :attr:`VDWFunc.param_names`, and combined_params is a list of
        combined VDW pair param values corresponding (in sequential order) to
        the properties :attr:`VDWFunc.pair_param_names`. The registry is
        populated by default with 'arithmetic/geometric', 'geometric', and
        'lb/geometric'.

        Returns: { str: f([float,...,float],[float,...,float],float) -> [float,...,float], ... }
        """
        return cls._comb_rule_map

    @classmethod
    def AddVDWCombRule(cls, key, comb_rule):
        """Add a new VDW combine rule to the registry.

        Arguments:
            key -- str

            comb_rule -- f([float,...,float],[float,...,float],float) -> [float,...,float]

        """
        cls._comb_rule_map[key] = comb_rule
        registry = _viparr.Rules.VDWCombRuleRegistry
        registry[key] = _viparr.Rules.VDWCombRule(comb_rule)
        _viparr.Rules.VDWCombRuleRegistry = registry

    @classmethod
    def DelVDWCombRule(cls, key):
        """Delete a VDW combine rule from the registry.

        Arguments:
            key -- str

        """
        del cls._comb_rule_map[key]
        registry = _viparr.Rules.VDWCombRuleRegistry
        del registry[key]
        _viparr.Rules.VDWCombRuleRegistry = registry

    def __init__(self):
        """Construct an empty :class:`Rules` object.
        
        Exclusion rule defaults to 1 and boolean flags take default values.
        """
        self._Rules = _viparr.Rules()

    def __eq__(self, other):
        try:
            return self._Rules == other._Rules
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._Rules.__hash__()

    @property
    def info(self):
        """List of info strings."""
        return self._Rules.info

    @info.setter
    def info(self, val):
        self._Rules.info = val

    @property
    def vdw_func(self):
        """VDW functional form.
        
        Should correspond to a key in :func:`Rules.GetVDWFuncRegistry()`, or be
        '' if this forcefield is to take the functional form of another
        forcefield during viparr execution.
        """
        return self._Rules.vdw_func

    @vdw_func.setter
    def vdw_func(self, val):
        self._Rules.vdw_func = val

    @property
    def vdw_comb_rule(self):
        """VDW combine rule.

        Should correspond to a key in :func:`Rules.GetVDWCombRuleRegistry()`,
        or be '' if this forcefield is to take the combine rule of another
        forcefield during viparr execution.
        """
        return self._Rules.vdw_comb_rule

    @vdw_comb_rule.setter
    def vdw_comb_rule(self, val):
        self._Rules.vdw_comb_rule = val

    @property
    def plugins(self):
        """List of plugins names.

        Each plugin name should correspond to a key in
        :func:`Forcefield.GetPluginRegistry()`.
        """
        return self._Rules.plugins

    @plugins.setter
    def plugins(self, val):
        self._Rules.plugins = val

    @property
    def fatal(self):
        """Whether viparr should fail if a required parameter is missing.
        
        Default is True.
        """
        return self._Rules.fatal

    @fatal.setter
    def fatal(self, val):
        self._Rules.fatal = val

    @property
    def nbfix_identifier(self):
        """A string identifier for the NBFix group of this forcefield.

        'nonbonded' parameters from all forcefields with the same
        nbfix_identifier are processed together during NB-fix.
        """
        return self._Rules.nbfix_identifier

    @nbfix_identifier.setter
    def nbfix_identifier(self, val):
        self._Rules.nbfix_identifier = val

    def setExclusions(self, exclusions, es_scale, lj_scale):
        """Set the exclusion rule and ES and LJ scale factors.

        The 'exclusions' argument indicates the max exclusion rule and should be
        2, 3, or 4, corresponding to 1-2, 1-3, or 1-4 exclusions. An
        'exclusions' argument of 1 indicates no exclusions except for those
        specified in the templates or by SmartsExclusions. 'es_scale' and
        'lj_scale' are of the form [scale for 1-2, scale for 1-3, ...] and
        should be of length 'exclusions'-1 to provide a scaling factor for each
        exclusion up to and including the max exclusion rule.

        Arguments:
            exclusions -- int

            es_scale -- [float, ..., float]

            lj_scale -- [float, ..., float]

        """
        self._Rules.setExclusions(exclusions, es_scale, lj_scale)

    @property
    def exclusions(self):
        """The max exclusion rule.
        
        A value of n indicates exclusions up to 1-n separation.
        """
        return self._Rules.exclusions

    def es_scale(self, separation):
        """Get the ES scaling factor for a given separation.

        Arguments:
            separation -- int (2, 3, or 4)

        Returns: float
        """
        return self._Rules.es_scale(separation)

    def lj_scale(self, separation):
        """Get the LJ scaling factor for a given separation.

        Arguments:
            separation -- int (2, 3, or 4)

        Returns: float
        """
        return self._Rules.lj_scale(separation)

    def __str__(self):
        output = ''
        if len(self.info) > 0:
            output += 'Info:\n'
            for line in self.info:
                output += line + '\n'
            output += '\n'
        output += 'VDW functional form: ' + self.vdw_func + '\n'
        output += 'VDW combine rule: ' + self.vdw_comb_rule + '\n'
        if self.nbfix_identifier != '':
            output += 'NB-fix identifier: ' + self.nbfix_identifier + '\n'
        if not self.fatal:
            output += 'Fatal: False\n'
        output += '\n'
        if len(self.plugins) > 0:
            output += 'Plugins:\n'
            for p in self.plugins:
                output += p + '\n'
            output += '\n'
        return output

class TemplateTyper(object):
    """Performs template matching and atom-typing."""

    @classmethod
    def _from_boost(cls, _TemplateTyper):
        if _TemplateTyper is None:
            return None
        typer = cls()
        typer._Typer = _TemplateTyper
        return typer

    def __init__(self):
        """Construct an empty typer with no templates."""
        self._Typer = _viparr.TemplateTyper()

    def __eq__(self, other):
        try:
            return self._Typer == other._Typer
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._Typer.__hash__()

    def addTemplate(self, tpl):
        """Add a template.

        Arguments:
            tpl -- :class:`TemplatedSystem`

        """
        self._Typer.addTemplate(tpl._TemplatedSystem)

    def delTemplate(self, tpl):
        """Delete a template.

        Arguments:
            tpl -- :class:`TemplatedSystem`

        """
        self._Typer.delTemplate(tpl._TemplatedSystem)

    def findTemplate(self, name):
        """Find all templates with the given name.

        Arguments:
            name -- str

        Returns: [:class:`TemplatedSystem`, ..., :class:`TemplatedSystem`]
        """
        return [TemplatedSystem._from_boost(tpl) for tpl in
                self._Typer.findTemplate(name)]

    def matchFragment(self, tsystem, atoms):
        """Match templates to a fragment of a system.

        'atoms' must correspond to a single complete fragment of 'tsystem'
        (which can be obtained by :func:`tsystem.system.updateFragids()`).
        Templates are matched separately to the separate residues within 
        'atoms'; matches are based on isomorphism of the bond graph.

        If a matching template is found for all residues in 'atoms', then the
        output is (matches, ""), where matches is a list of (tpl, atom_list)
        pairs corresponding to the matched templates and atoms in 'tsystem'
        matched by the templates. If a residue in 'atoms' cannot be matched,
        then the output is ([], err_msg) where err_msg indicates a possible
        reason why the matching failed.

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: ([(:class:`TemplatedSystem`, [:class:`msys.Atom`, ... ]), ... ], str)
        """
        out, msg = self._Typer.matchFragment(tsystem._TemplatedSystem,
                [a.id for a in atoms])
        matches = [(TemplatedSystem._from_boost(pair[0]),
            [msys.Atom(tsystem.system._ptr,
                id) for id in pair[1]]) for pair in out]
        return (matches, msg)

    def assignMatch(self, tsystem, matches, rename_atoms=False,
            rename_residues=False):
        """Map information from a list of matched templates to a system.

        'matches' should be the output from :func:`matchFragment`. Permanent
        changes to :attr:`tsystem.system` are the mapping of atom charges from
        the template, the addition of pseudo particles and pseudo bonds, and
        optionally the mapping of atom names and residue names. Other changes to
        tsystem (but not to the contained :attr:`tsystem.system`) are the
        addition of atom type and pset information, the possible addition of
        bond aromaticity, and the addition of lists of typed atoms, non-pseudo
        bonds, pseudo-bonds, angles, dihedrals, impropers, extra exclusions,
        cmaps, and pseudo-sites to be matched during parameter matching.
        
        Arguments:
            tsystem -- :class:`TemplatedSystem`

            matches -- [(:class:`TemplatedSystem`, [:class:`msys.Atom`, ... ]), ... ]

            rename_atoms -- bool

            rename_residues -- bool

        """
        input = [(pair[0]._TemplatedSystem,
            [a.id for a in pair[1]]) for pair in matches]
        self._Typer.assignMatch(tsystem._TemplatedSystem, input, rename_atoms,
                rename_residues)

    @property
    def templates(self):
        """List of all contained templates."""
        return [TemplatedSystem._from_boost(tpl) \
                for tpl in self._Typer.templates()]

    def __str__(self):
        output = 'Templates:\n'
        for template in self.templates:
            output += template.system.residue(0).name + ': '
            output += template.__repr__() + '\n'
        output += '\n'
        return output

###################### Import, export, and merge ###############################
def ImportForcefield(dir, require_rules=True):
    """Import entire forcefield from a forcefield directory.
    
    :class:`Rules` is imported from 'rules' and templates from any files of
    name 'templates...'. If 'atomtypes.def' exists, then a :class:`SmartsTree`
    or :class:`ScoredSmarts` object is imported from 'atomtypes.def' and 
    (optional) 'extendtypes.def', :class:`SmartsExclusions` are imported from 
    (optional) 'exclusions.def', :class:`SmartsImpropers` are imported from 
    (optional) 'impropers.def', :class:`msys.ParamTables` are imported from 
    (optional) 'virtuals.def' and 'drudes.def', and the :class:`Forcefield` is 
    constructed with a :class:`SmartsTyper` or :class:`ScoredSmartsTyper` 
    (depending on the 'score_based_atomtyping' flag in rules).  Otherwise, the
    :class:`Forcefield` is constructed with a :class:`TemplateTyper`. If
    require_rules is True, an error is raised if the rules file is missing. 

    'cmap' is treated as a cmap param file, and all other files are treated as
    regular param files. Each param file is imported into the static dictionary
    of shared param tables in the :class:`Forcefield` class, with the list of
    added parameters stored with the :class:`Forcefield` object. The list of
    parameters is always stored with the forcefield object under the same name
    as the param file; the static shared param table to which the parameters are
    imported is the table with key being the param file name.

    Arguments:
        dir -- str

        require_rules -- bool

    Returns: :class:`Forcefield`

    Side effect: Modifies the static param table dictionary of the
        :class:`Forcefield` class

    """
    return Forcefield._from_boost(_viparr.ImportForcefield(dir, require_rules,
        {}))

def ImportRules(path):
    """Import a rules file.
    
    All fields of the rules file are optional. If 'exclusions' is unspecified,
    max exclusion rule is set to length of 'es_scale' + 1 (or length of
    'lj_scale' + 1; these must be equal), or to a default value of 4 if neither
    scale is specified. Scales default to 0 for all interaction distances if
    unspecified. Boolean flags take their default values if unspecified; see
    documentation on the Rules class for details.

    Arguments:
        path -- str

    Returns: :class:`Rules`
    """
    return Rules._from_boost(_viparr.ImportRules(path))

def ImportTemplates(path):
    """Import a template file.
    
    Bonds, impropers, etc. can reference external atom names (usually indicated
    by $num), but the reference names must be consistent. Atom types cannot
    contain whitespace characters.
    
    Arguments:
        path -- str

    Returns: [:class:`TemplatedSystem`, ..., :class:`TemplatedSystem`]
    """
    return [TemplatedSystem._from_boost(tpl) \
            for tpl in _viparr.ImportTemplates(path)]

def ImportCmap(path):
    """Import a cmap parameter file.

    Arguments:
        path -- str

    Returns: msys.ParamTable
    """
    return [msys.ParamTable(table) for table in _viparr.ImportCmap(path)]

def ImportParams(table_name, path, nbfix_identifier=''):
    """Import a parameter file.
    
    If 'table_name' exists as a key in the :class:`Forcefield` class' static
    dictionary of param tables, then new params are appended to this table.
    Otherwise, a new table is created and added to the dictionary with key
    'table_name'. Additional "type" and "memo" columns, as well as a
    "nbfix_identifier" column if 'table_name' is vdw1 or vdw2, are added. The
    "type" column is a ' '-concatenation of all of the "type" fields of the
    file. The "nbfix_identifier" column for vdw1 and vdw2 tables is populated
    with the given value. A list of added params is returned.

    "type" fields cannot have whitespace characters. A "params" key of "cmap" is
    treated specially as a cmap ID; a value of 1 in the "cmap" column is stored
    as a string "cmap1" in the table, etc. "nbfix_identifier" cannot be empty
    for vdw2 tables.

    Arguments:
        table_name -- str

        path -- str

        nbfix_identifier -- str

    Returns: [:class:`msys.Param`, ..., :class:`msys.Param`]

    Side effect: Modifies the static dictionary of param tables in the
        :class:`Forcefield` class

    """
    rowIDs = _viparr.ImportParams(table_name, path, [], nbfix_identifier)
    table = Forcefield.ParamTable(table_name)
    return [msys.Param(table._ptr, id) for id in rowIDs]

def ExportForcefield(ff, dir):
    """Export entire forcefield to an empty directory.
    
    :class:`Rules` are exported to 'rules', templates to 'templates', 
    :class:`SmartsTree` to 'atomtypes.def', :class:`SmartsExclusions` to
    'exclusions.def', :class:`SmartsImpropers` to 'impropers.def', virtuals and
    drudes tables to 'virtuals.def' and 'drudes.def', cmap tables to 'cmap', and
    the remaining param tables to the table name in the forcefield. Templates
    are not exported with a :class:`SmartsTyper`.

    Arguments:
        ff -- Forcefield

        dir -- str, path to non-existent or empty directory

    """
    _viparr.ExportForcefield(ff._Forcefield, dir)

def ExportRules(rules, path):
    """Export a :class:`Rules` object to a file.

    Arguments:
        rules -- :class:`Rules`

        path -- str

    """
    _viparr.ExportRules(rules._Rules, path)

def ExportTemplates(templates, path):
    """Export all templates to a file.

    All templates are exported to a single file, even if they were imported from
    multiple files.

    Arguments:
        templates -- [:class:`TemplatedSystem`, ..., :class:`TemplatedSystem`]

        path -- str

    """
    _viparr.ExportTemplates([tpl._TemplatedSystem for tpl in templates], path)

def ExportCmap(cmap_tables, path):
    """Export cmap parameter tables to a file.

    Arguments:
        cmap_tables -- [:class:`msys.ParamTable`, ..., :class:`msys.ParamTable`]

        path -- str

    """
    _viparr.ExportCmap([table._ptr for table in cmap_tables], path)

def ExportParams(params, path):
    """Export a subset of parameters in a parameter table to a file.

    'params' must belong to the same param table.

    Arguments:

        params -- [:class:`msys.Param`, ..., :class:`msys.Param`]

        path -- str

    """
    if len(params) == 0:
        return
    table = params[0]._ptr
    _viparr.ExportParams(table, [param.id for param in params], path)

def MergeForcefields(src, patch, append_only=False, verbose=True):
    """Merge forcefield patch into forcefield src.

    Merges rules, templates, param tables, scored SMARTS atomtypes, SMARTS
    exclusions, SMARTS impropers, virtual definition tables, and drude
    definition tables. Merging is not supported for SMARTS atomtype trees and
    cmap tables; if these are present in patch, they directly overwrite those in
    src. See individual merge functions for additional details.

    Will throw an exception if append_only is True and the merge procedure
    attempts to overwrite any components of src.

    Arguments:
        src -- :class:`Forcefield`

        patch -- :class:`Forcefield`

        append_only -- bool

        verbose -- bool

    """
    _viparr.MergeForcefields(src._Forcefield, patch._Forcefield, append_only,
            verbose)

def MergeRules(src_rules, patch_rules, verbose=True):
    """Merge rules and plugins of patch into src.

    VDW functional form and combine rule, exclusion rule, and ES and LJ scaling
    factors, if defined in patch, must agree with those in src.

    Arguments:
        src_rules -- :class:`Rules`

        patch_rules -- :class:`Rules`

        verbose -- bool

    """
    _viparr.MergeRules(src_rules._Rules, patch_rules._Rules, verbose)

def MergeTemplates(src_typer, patch_typer, append_only=False, verbose=True):
    """Merge templates of patch into src.

    If append_only is False, templates of patch overwrite any templates in src
    of the same name.

    Arguments:
        src_typer -- :class:`TemplateTyper`

        patch_typer -- :class:`TemplateTyper`

        append_only -- bool

        verbose -- bool

    """
    _viparr.MergeTemplates(src_typer._Typer, patch_typer._Typer, append_only,
            verbose)

def MergeParams(src_params, patch_params, append_only=False, verbose=True):
    """Merge parameters of patch into src for a single param table.

    If append_only is False, rows in patch overwrite any rows in src matching
    the same "type" field. In the merged table, rows in patch not present in
    src appear first, in their original order, followed by rows in src, in
    their original order and possibly with parameter values overwritten by
    those in patch.

    Arguments:
        src_params -- [:class:`msys.Param`, ..., :class:`msys.Param`]

        patch_params -- [:class:`msys.Param`, ..., :class:`msys.Param`]

        append_only -- bool

        verbose -- bool

    """
    if len(patch_params) == 0:
        return src_params
    ids = _viparr.MergeParams([p.id for p in src_params],
            [p.id for p in patch_params], patch_params[0].table._ptr,
            append_only, verbose)
    return [patch_params[0].table.param(id) for id in ids]

######################## Parameter matching classes ############################
class Pattern(object):
    """Structure that holds atom types, bond types, and string flags to match.

    :class:`Pattern` is a helper class for :class:`ParameterMatcher`; see 
    :class:`ParameterMatcher` for additional documentation. Flags are currently
    only used for pseudo particles, to match 'pset'. Bond type conventions:

        '-' : order-1 bond

        '=' : order-2 bond

        '#' : order-3 bond

        ':' : aromatic bond

    """

    @classmethod
    def _from_boost(cls, _Pattern):
        if _Pattern is None:
            return None
        pattern = cls()
        pattern._Pattern = _Pattern
        return pattern

    def __init__(self):
        """Construct empty pattern with no atoms, bonds, or flags."""
        self._Pattern = _viparr.Pattern()

    def __eq__(self, other):
        try:
            return self._Pattern == other._Pattern
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._Pattern.__hash__()

    @property
    def atoms(self):
        """List of atom types to match."""
        return self._Pattern.atoms

    @atoms.setter
    def atoms(self, val):
        self._Pattern.atoms = val

    @property
    def bonds(self):
        """List of bond types to match."""
        return self._Pattern.bonds

    @bonds.setter
    def bonds(self, val):
        self._Pattern.bonds = val

    @property
    def flags(self):
        """List of additional string flags to match."""
        return self._Pattern.flags

    @flags.setter
    def flags(self, val):
        self._Pattern.flags = val

    def __str__(self):
        return self._Pattern.__str__()

class SystemToPattern(object):
    """Functions to create a :class:`Pattern` from an atom-tuple."""

    @staticmethod
    def NBType(tsystem, atoms):
        """Create a :class:`Pattern` using the atoms' non-bonded type.

        'atoms' should belong to tsystem. The returned :class:`Pattern` is:

            atom types: nbtype(a_0), ..., nbtype(a_n)

            bond types: None

            flags: None

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.SystemToPattern.NBType(
            tsystem._TemplatedSystem, [atom.id for atom in atoms]))

    @staticmethod
    def BType(tsystem, atoms):
        """Create a :class:`Pattern` using the atoms' bonded type.

        'atoms' should belong to tsystem. The returned :class:`Pattern` is:

            atom types: btype(a_0), ..., btype(a_n)

            bond types: None

            flags: None

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.SystemToPattern.BType(
            tsystem._TemplatedSystem, [atom.id for atom in atoms]))

    @staticmethod
    def Bonded(tsystem, atoms):
        """Create a :class:`Pattern` for atoms in a bonded chain.

        'atoms' should belong to tsystem and be bonded in their given order.
        The returned :class:`Pattern` is:

            atom types: btype(a_0), ..., btype(a_n)

            bond types: bond_order(a_0,a_1), ..., bond_order(a_{n-1},a_n)

            flags: None

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.SystemToPattern.Bonded(
            tsystem._TemplatedSystem, [atom.id for atom in atoms]))

    @staticmethod
    def BondToFirst(tsystem, atoms):
        """Create a :class:`Pattern` for atoms bonded to the first atom.

        'atoms' should belong to tsystem and be bonded to the first atom.
        The returned :class:`Pattern` is:

            atom types: btype(a_0), ..., btype(a_n)

            bond types: bond_order(a_0,a_1), bond_order(a_0,a_2), ..., bond_order(a_0,a_n)

            flags: None

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.SystemToPattern.BondToFirst(
            tsystem._TemplatedSystem, [atom.id for atom in atoms]))

    @staticmethod
    def PseudoBType(tsystem, atoms):
        """Create a :class:`Pattern` from pseudo 'pset' and types of site atoms.

        'atoms' should belong to tsystem, where the first atom is a pseudo
        and the remaining atoms are its other site atoms.
        The returned :class:`Pattern` is:

            atom types: btype(a_1), ..., btype(a_n)

            bond types: None

            flags: pset(a_0)

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.SystemToPattern.PseudoBType(
            tsystem._TemplatedSystem, [atom.id for atom in atoms]))

    @staticmethod
    def PseudoBondToFirst(tsystem, atoms):
        """Create a :class:`Pattern` from pseudo 'pset' and bonds to first site.

        'atoms' should belong to tsystem, where the first atom is a pseudo, the
        remaining atoms are its other site atoms, and all other site atoms are
        bonded to the first site atom (usually the parent of the pseudo).
        The returned :class:`Pattern` is:

            atom types: btype(a_1), ..., btype(a_n)

            bond types: bond_order(a_1,a_2), bond_order(a_1,a_3), ..., bond_order(a_1,a_n)

            flags: pset(a_0)

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.SystemToPattern.PseudoBondToFirst(
            tsystem._TemplatedSystem, [atom.id for atom in atoms]))

class TypeToPattern(object):
    """Functions to create a :class:`Pattern` from 'type' column of param table.
    """

    @staticmethod
    def Default(type_string):
        """Create a :class:`Pattern` from atom types and bond orders.

        Tokenizes 'type_string' using ' ', treats tokens '-', '=', '#', and ':',
        and '~' as bond types, and treats all remaining tokens as atom types.

        Arguments:
            type_string -- str

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.TypeToPattern.Default(type_string))

    @staticmethod
    def Pseudo(type_string):
        """Create a :class:`Pattern` from atom types, bond orders, and pset.

        Tokenizes 'type_string' using ' ', treats tokens '-', '=', '#', and ':',
        and '~' as bond types, treats the last token (assumed to be pset) as
        a flag, and treats all remaining tokens as atom types.

        Arguments:
            type_string -- str

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.TypeToPattern.Pseudo(type_string))

class Permutation(object):
    """Functions that return a permutation of a :class:`Pattern` object."""

    @staticmethod
    def Identity(pattern):
        """Does not modify 'pattern'.
        
        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Identity(
            pattern._Pattern))

    @staticmethod
    def Reverse(pattern):
        """Reverses the order of atoms and the order of bonds in 'pattern'.

        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Reverse(
            pattern._Pattern))

    @staticmethod
    def Improper1(pattern):
        """One of six functions that generates permutations for impropers.

        Assumes that 'pattern' has four atoms and three bonds, fixes the first
        atom, and permutes the remaining three atoms and corresponding bonds.
        :func:`Identity`, :func:`Improper1`, :func:`Improper2`,
        :func:`Improper3`, :func:`Improper4`, and :func:`Improper5`
        provide the six different permutations of these three atoms/bonds.

        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Improper[1](
            pattern._Pattern))

    @staticmethod
    def Improper2(pattern):
        """One of six functions that generates permutations for impropers.

        See :func:`Improper1` documentation.

        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Improper[2](
            pattern._Pattern))

    @staticmethod
    def Improper3(pattern):
        """One of six functions that generates permutations for impropers.

        See :class:`Improper1` documentation.

        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Improper[3](
            pattern._Pattern))

    @staticmethod
    def Improper4(pattern):
        """One of six functions that generates permutations for impropers.

        See :class:`Improper1` documentation.

        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Improper[4](
            pattern._Pattern))

    @staticmethod
    def Improper5(pattern):
        """One of six functions that generates permutations for impropers.

        See :class:`Improper1` documentation.

        Arguments:
            pattern -- :class:`Pattern`

        Returns: :class:`Pattern`
        """
        return Pattern._from_boost(_viparr.Permutation.Improper[5](
            pattern._Pattern))

class ParameterMatcher(object):
    """Matches atom tuples to a subset of rows in a param table.
    
    A :class:`ParameterMatcher` stores a set of parameters from a param
    table of a forcefield and can match atom tuples from a system against the
    param table's 'type' column. Almost all of viparr's built-in forcefield
    plugins use this :class:`ParameterMatcher` class to perform parameter
    matching; flexibility is provided in the constructor to control how tuples
    and rows are to be matched.
    
    This class supports matching using hierarchical atom types, if atom types
    were generated using a :class:`SmartsTyper`. For atom types generated by a 
    :class:`TemplateTyper`, non-hierarchical matching is used.

    Typical use:
    matcher = ParameterMatcher.FromFF(ff, table_name, ...);
    param = matcher.match(tsystem, tuple);
    term_table.addTerm(tuple, param);
    """
    @classmethod
    def _from_boost(cls, _ParameterMatcher):
        if _ParameterMatcher is None:
            return None
        matcher = cls([], SystemToPattern.BType, TypeToPattern.Default,
                [Permutation.Identity])
        matcher._ParameterMatcher = _ParameterMatcher
        return matcher

    # We do not expose the SystemToPattern, TypeToPattern, and Permutation
    # function-wrapper classes, but rather store a map from Python functions to
    # the wrapped functions and provide an interface to operate directly on
    # Python functions
    _stp_map = { SystemToPattern.NBType: _viparr.SystemToPattern.NBType, 
            SystemToPattern.BType: _viparr.SystemToPattern.BType,
            SystemToPattern.Bonded: _viparr.SystemToPattern.Bonded,
            SystemToPattern.BondToFirst: _viparr.SystemToPattern.BondToFirst,
            SystemToPattern.PseudoBType: _viparr.SystemToPattern.PseudoBType,
            SystemToPattern.PseudoBondToFirst:
                _viparr.SystemToPattern.PseudoBondToFirst }


    _ttp_map = { TypeToPattern.Default: _viparr.TypeToPattern.Default,
            TypeToPattern.Pseudo: _viparr.TypeToPattern.Pseudo }

    _perm_map = { Permutation.Identity: _viparr.Permutation.Identity,
            Permutation.Reverse: _viparr.Permutation.Reverse }

    def __init__(self, params, sys_to_pattern, type_to_pattern, permutations):
        """Construct from list of parameters, atom-types, and helper functions.

        The parameter matcher operates by converting both the atom-tuple and
        the 'type' field of the parameter into :class:`Pattern` objects and
        matching these :class:`Pattern` objects. The functions that convert an
        atom-tuple to a :class:`Pattern` and a 'type' field to a pattern are
        user-supplied; some pre-defined functions that do these conversions are
        found in the :class:`SystemToPattern` and :class:`TypeToPattern`
        classes. The matcher can also match various permutations of the
        atom-tuple :class:`Pattern` object to the 'type' :class:`Pattern`
        object; these permutation functions are again user-supplied, with some
        pre-defined permuations available in the :class:`Permutation` class.

        Arguments:
            params -- [:class:`msys.Param`, ..., :class:`msys.Param`]

            sys_to_pattern -- f([:class:`msys.Atom`, ... ]) -> :class:`Pattern`

            type_to_pattern -- f(str) -> :class:`Pattern`

            permutations -- [f(:class:`Pattern`) -> :class:`Pattern`, ... ]

        """
        if sys_to_pattern in self.__class__._stp_map:
            stp = self.__class__._stp_map[sys_to_pattern]
        else:
            def f(_TemplatedSystem, ids):
                return sys_to_pattern(
                        TemplatedSystem._from_boost(_TemplatedSystem),
                        [msys.Atom(_TemplatedSystem.system(),
                            id) for id in ids])._Pattern
            stp = _viparr.SystemToPattern(f)
            self.__class__._stp_map[sys_to_pattern] = stp
        if type_to_pattern in self.__class__._ttp_map:
            ttp = self.__class__._ttp_map[type_to_pattern]
        else:
            def f(type_string):
                return type_to_pattern(type_string)._Pattern
            ttp = _viparr.TypeToPattern(f)
            self.__class__._ttp_map[type_to_pattern] = ttp
        perms = []
        for permutation in permutations:
            if permutation in self.__class__._perm_map:
                perms.append(self.__class__._perm_map[permutation])
            else:
                def f(_Pattern):
                    return permutation(Pattern._from_boost(_Pattern))._Pattern
                perm = _viparr.Permutation(f)
                self.__class__._perm_map[permutation] = perm
                perms.append(perm)
        if len(params) == 0:
            table = msys.CreateParamTable()._ptr
        else:
            table = params[0]._ptr
        self._ParameterMatcher = _viparr.ParameterMatcher(table,
                [param.id for param in params], stp, ttp, perms)

    @classmethod
    def FromFF(cls, forcefield, table_name, sys_to_pattern, type_to_pattern,
            permutations):
        """Construct from a :class:`Forcefield` and table name.

        Calls the constructor using :func:`forcefield.params(table_name)` as
        the list of parameters. See the constructor for documentation on the
        other arguments.

        Arguments:
            forcefield -- :class:`Forcefield`

            table_name -- str

            sys_to_pattern -- f([:class:`msys.Atom`, ... ]) -> :class:`Pattern`

            type_to_pattern -- f(str) -> :class:`Pattern`

            permutations -- [f(:class:`Pattern`) -> None, ... ]

        """
        params = forcefield.params(table_name)
        matcher = cls(params, sys_to_pattern, type_to_pattern, permutations)
        return matcher

    def __eq__(self, other):
        try:
            return self._ParameterMatcher == other._ParameterMatcher
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._ParameterMatcher.__hash__()

    def match(self, tsystem, atoms, allow_repeat=False):
        """Matches a single atom tuple to the contained parameters.

        Constructs an atom-tuple :class:`Pattern` object and finds a matching
        'type' :class:`Pattern` object. Returns the matching parameter and
        permutation function if a match is found, or (invalid_param, None)
        otherwise, where 'invalid_param' has ID -1. If allow_repeat is True,
        multiple matches of identical types are allowed; otherwise an error is
        thrown.

        A 'type' :class:`Pattern` object 'tpattern' is a match for the 
        atom-tuple :class:`Pattern` object 'apattern' if there exists a valid 
        permutation 'apatternperm' of 'apattern' such that all of the following 
        are true:

        1. tpattern.atoms equals apatternperm.atoms, or (in the case of
            hierarchical matching) each atom type in tpattern.atoms has a
            descendant in the type-hierarchy tree that equals the corresponding
            atom type in apatternperm.atoms.

        2. tpattern.bonds is empty, or tpattern.bonds equals apatternperm.bonds

        3. tpattern.flags equals apatternperm.flags

        Note that if the 'type' :class:`Pattern` constructed from the forcefield
        parameter table does not specify bond types, then bonds are ignored even
        if they are present in the atom-tuple :class:`Pattern`. '*' is used as a
        wild-character for atom types and '~' is used as a wild-character for
        bond types; 'equals' in the conditions above means with respect to these
        wild-characters (e.g. atom-type 'c*' equals atom-type 'car').

        In the case of multiple matching 'type' patterns, a single best match is
        chosen by the following rules:

        1. A match without a '*' wild-character in the atom type is better
           than a match with a '*' wild-character.

        2. Among matches without a '*' wild-character, a match whose atom types
            have higher hierarchical priority is better than a match whose atom
            types have lower hierarchical priority. (For non-hierarchical
            matching, all non-'*' matches are at the same priority.) If more
            than one 'type' Pattern matches at the same priority, an error is
            raised unless the matched 'type' patterns are identical and
            allow_repeat is True.

        3. Among matches with a '*' wild-character, the first match in the
            table is best, regardless of priority. No error is raised for
            multiple matches.

        In hierarchical matching, the priority for a tuple of atom types is a
        tuple of integers, and it is defined as follows:
        (0,0,...,0) - all atom types are root types
        (0,0,...,1) - exactly one atomtype is one level below its root
        (0,0,...,2) - exactly one atomtype is two levels below its root
        (0,...,0,1,1) - exactly two atomtypes are one level below their roots 
        etc.

        Priorities are ordered in standard lexicographical order, with
        (0,0,...,0) being the lowest priority. Hence the highest priority match
        is the match that maximizes the closest distance of any atom type in
        the matched tuple to its root type in the tree.

        Arguments:
            tsystem -- :class:`TemplatedSystem`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        Returns: (:class:`msys.Param`, f(:class:`Pattern`) -> None)
        """
        row_id, perm = self._ParameterMatcher.match(
                tsystem._TemplatedSystem, [atom.id for atom in atoms],
                allow_repeat)
        for py_perm, boost_perm in self.__class__._perm_map.items():
            if boost_perm == perm:
                return (msys.Param(self._ParameterMatcher.paramTable(), row_id),
                        py_perm)
        if perm == _viparr.Permutation.Null:
            return (msys.Param(self._ParameterMatcher.paramTable(), -1), None)
        raise RuntimeError('VIPARR bug -- permutation not found')

    def writeMultiple(self, param, atoms, term_table):
        """Writes matches to a system for param tables with duplicate rows.

        Only used in the 'propers' plugin to handle some forcefields where
        multiple rows of dihedral_trig parameters are defined for a single type
        pattern.

        For a given tuple of atoms and matched parameter, all parameters with
        identical 'type' value as the matched parameter are found, and one entry
        is written to the given term table for each such match.

        Arguments:
            param -- :class:`msys.Param`

            atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

            term_table -- :class:`msys.TermTable`

        """
        self._ParameterMatcher.writeMultiple(param.id,
                [atom.id for atom in atoms], term_table._ptr)

    @property
    def params(self):
        """The parameters being matched."""
        table = self._ParameterMatcher.paramTable()
        return [msys.Param(table, id) for id in self._ParameterMatcher.rowIDs()]

    @property
    def sys_to_pattern(self):
        """The system-to-pattern function used to construct the matcher."""
        stp = self._ParameterMatcher.sysToPattern()
        for key, item in self.__class__._stp_map.items():
            if item == stp:
                return key
        raise RuntimeError('VIPARR error -- sys_to_pattern not found')

    @property
    def type_to_pattern(self):
        """The type-to-pattern function used to construct the matcher."""
        ttp = self._ParameterMatcher.typeToPattern()
        for key, item in self.__class__._ttp_map.items():
            if item == ttp:
                return key
        raise RuntimeError('VIPARR bug -- type_to_pattern not found')

######################### High-level functions #################################
def AddSystemTables(system, table_name_mapping=None):
    """Merge parameters from a system into the static shared param tables.

    For each term table in the system, adds the parameters of the associated
    param table to the corresponding table in the static dictionary of shared
    param tables of the :class:`Forcefield` class, and associates the term table
    with the shared param table. 

    In most cases, parameters for a term table named 'some_table' are added to
    the shared param table :attr:`Forcefield.ParamTable('some_table')` of the
    same name.  If the term table name contains any key of the given
    'table_name_mapping' dictionary as a substring, then parameters are added to
    the shared param table with a modified name, where the first instance of the
    contained key is replaced with its mapped value in 'table_name_mapping'. If
    'table_name_mapping' is None, the default mapping used is
    {'nonbonded': 'vdw1', 'virtual': 'virtuals'}

    Arguments:
        system -- :class:`msys.System`

        table_name_mapping -- {str: str, ..., str: str}

    Side effect: Modifies the static dictionary returned by
        :func:`Forcefield.AllParamTables()`

    """
    _viparr.AddSystemTables(system._ptr, {}, table_name_mapping)

def ExecuteViparr(system, ffs, atoms=None, rename_atoms=False,
        rename_residues=False, with_constraints=True, fix_masses=True,
        fatal=True, compile_plugins=True, verbose=True, verbose_matching=True):
    """Run viparr to parametrize a system using a list of forcefields.

    Equivalent to the viparr command-line executable without reorder-ids
    or forcefield merging. (Forcefields may be merged first using
    :func:`MergeForcefields`.) Each fragment of the system is parametrized
    using the first forcefield that can match all atoms/templates in
    the fragment.  Options are provided to parametrize only a subset of
    atoms, to copy atom and residue names from the template definitions
    for template-based forcefields, and to build constraints and fix
    atom masses at the end of viparr execution.  Charges, masses, BCI
    charge corrections, scaled pairs terms, and virtuals_shift parameters
    are added to the system, and temporary tables added by viparr are
    removed unless 'compile_plugins' is set to False. If 'fatal' is False,
    the :attr:`Rules.fatal` flag of all forcefields are set to False.

    If the system has existing parameter tables, they will be merged
    into the static shared parameter tables of the :class:`Forcefield`
    class using :func:`AddSystemTables`, and parameters for 'atoms'
    will be overwritten.  'atoms' must correspond to a set of complete
    bonded fragments.

    Original atom IDs are preserved and IDs of added pseudo particles are at
    the end; to reorder the IDs so that pseudos are next to parents, run
    :func:`ReorderIDs`.

    Arguments:
        system -- :class:`msys.System`

        ffs -- [:class:`Forcefield`, ..., :class:`Forcefield`]

        atoms -- [:class:`msys.Atom`, ... ], default: :attr:`system.atoms`

        rename_atoms -- bool

        rename_residues -- bool

        with_constraints -- bool

        fix_masses -- bool

        fatal -- bool

        compile_plugins -- bool

        verbose -- bool

    """
    if atoms is None:
        atoms = system.atoms
    _viparr.ExecuteViparr(system._ptr, [ff._Forcefield for ff in ffs],
            [atom.id for atom in atoms], rename_atoms, rename_residues,
            with_constraints, fix_masses, fatal, compile_plugins, verbose, verbose_matching)

class CompilePlugins(object):
    """A collection of plugin compilation functions and helper functions.
    
    All functions in this class, with the exception of CleanupSystem, can be
    applied multiple times to a system to update its derived properties and
    parameters after parameters in the raw parameter tables change. Once
    CleanupSystem is performed, the unnecessary raw parameter tables are
    removed and plugin compilation can no longer be performed.
    """

    @staticmethod
    def CompilePlugins(system, plugins=None):
        """Compile all default plugins, or the list of plugins provided.

        By default, compiles masses, charges, BCI, improper trig,
        ureybradley, exchange-dispersion, pairs, and virtuals shift.

        Arguments:
            system -- :class:`msys.System`
            plugins -- [string, ...], eg. ["mass", "charges_formal"] or ff.rules.plugins where ff is a viparr.Forcefield object
        """
        default_list = ["mass", "charges_formal", "impropers", 
                        "ureybradley", "exclusions_and_scaled_pairs", 
                        "pairs_lj_scaled_14"]

        if plugins is None:
            plugins = default_list
                    
        _viparr.CompilePlugins(system._ptr, plugins)

    @staticmethod
    def CompileMasses(system):
        """Copy values from 'mass' table to atom mass property.

        Arguments:
            system -- :class:`msys.System`

        """
        #_viparr.CompileMasses(system._ptr)
        Forcefield.GetPluginRegistry()["mass"].compile(system)

    @staticmethod
    def CompileImproperTrig(system):
        """Copy 'improper_trig' terms and parameters to 'dihedral_trig' table.

        Arguments:
            system -- :class:`msys.System`

        """
        Forcefield.GetPluginRegistry()["impropers"].compile(system)

    @staticmethod
    def CompileUreyBradley(system):
        """Copy 'ureybradley_harm' terms and parameters to 'stretch_harm' table.

        Arguments:
            system -- :class:`msys.System`

        """
        Forcefield.GetPluginRegistry()["ureybradley"].compile(system)

    @staticmethod
    def ApplyNBFix(system):
        """Apply 'vdw2' term table params as overrides to 'nonbonded' table.

        Arguments:
            system -- :class:`msys.System`
        """
        _viparr.ApplyNBFix(system._ptr)

    @staticmethod
    def AddPairsTable(system):
        """Add a nonbonded-pairs term table to the system.

        Arguments:
            system -- :class:`msys.System`

        Returns: :class:`msys.TermTable`
        """
        return msys.TermTable(_viparr.AddPairsTable(system._ptr))

    @staticmethod
    def CompilePairsESScaled(system, pairs):
        """Add ES scaled pairs for non-zero scale factors.

        Arguments:
            system -- :class:`msys.System`
            pairs -- :class:`msys.TermTable`

        """
        Forcefield.GetPluginRegistry()["exclusions_and_scaled_pairs"].compile(system)

    @staticmethod
    def CompilePairsLJScaled(system, pairs):
        """Add LJ scaled pairs for non-zero scale factors.

        This automatically looks in the 'vdw1_14' table first, if present, for
        nonbonded parameters for 1-4 interactions.

        Arguments:
            system -- :class:`msys.System`
            pairs -- :class:`msys.TermTable`

        """
        Forcefield.GetPluginRegistry()["pairs_lj_scaled_14"].compile(system)

    @staticmethod
    def CompileScaledPairOverrides(system, pairs):
        """Moves 'scaled_pair_overrides' terms into the pairs table.

        Arguments:
            system -- :class:`msys.System`
            pairs -- :class:`msys.TermTable`

        """
        Forcefield.GetPluginRegistry()["scaled_pair_overrides"].compile(system)

    @staticmethod
    def CleanupSystem(system):
        """Remove all temporary tables added by viparr.

        Removes all term tables with no msys category or no terms.

        Arguments:
            system -- :class:`msys.System`

        """
        _viparr.CleanupSystem(system._ptr)

def BuildConstraints(system, atoms=None, keep=False, exclude=[], verbose=True):
    """Add constraints to the system.

    This is executed by default when viparr is evoked from the command line or
    when :func:`ExecuteViparr` is called.
    
    Adds constraint parameters to the system or a subset of atoms in the
    system, and updates "constrained" fields of stretch_harm and angle_harm
    tables. Supports HOH and AHn constraints. If keep is true, sets
    "constrained" to 0 for all stretch_harm and angle_harm params; otherwise,
    sets "constrained" to 1 for constrained params and 0 for unconstrained
    params. Optionally takes a set of constraint types to ignore; types can
    be "hoh", "ah1", "ah2", etc. 'system' must have has stretch_harm and
    angle_harm terms for all bonds and angles involving the specified atoms.
    If no atoms are specified, constraints are built for the entire system.

    Arguments:
        system -- :class:`msys.System`

        atoms -- [:class:`msys.Atom`, ..., :class:`msys.Atom`]

        keep -- bool

        exclude -- [str, ..., str]

        verbose -- bool

    """
    if atoms is None:
        atoms = system.atoms
    _viparr.BuildConstraints(system._ptr, [atom.id for atom in atoms], keep,
            exclude, verbose)

def FixMasses(system, atoms=None, verbose=True):
    """Equate masses for all atoms of the same element.

    This is executed by default when viparr is evoked from the command line or
    when :func:`ExecuteViparr` is called.
    
    Replaces the mass of all atoms of each element with the median mass of these
    atoms, to correct mass discrepancies between different forcefields. An
    option is provided to run over a subset of atoms.

    Arguments:
        system -- :class:`msys.System`

        atoms -- [:class:`msys.Atom`, ... ], default: :class:`system.atoms`

        verbose -- bool

    """
    if atoms is None:
        atoms = system.atoms
    _viparr.FixMasses(system._ptr, [atom.id for atom in atoms], verbose)

def RotateToMinAABB(system, cubic=True):
    """Find and apply the rigid rotation to the system that minimizes the
    volume of the bounding box around it. This box can be constrained
    to be cubic; otherwise, the height, width and length will
    generally be distinct numbers.

    Also translates the atoms in a system so that the centroid is at
    the origin.

    Arguments:
        system -- :class:`msys.System`
        cubic  -- bool

    """
    _viparr.RotateToMinAABB(system._ptr, cubic)

def ReorderIDs(system):
    """Reorder atom IDs in a parametrized system based on bond connectivity.

    Returns a clone of the original system with reordered IDs. Can be used after
    viparr parametrization to place pseudo IDs next to parents. Also reorders
    pair and exclusion terms so that the lower-ID atom always comes first.

    Arguments:
        system -- :class:`msys.System`

    Returns: :class:`msys.System`
    """
    return msys.System(_viparr.ReorderIDs(system._ptr))

def GenerateStandardPlugin(name, tuple_type, required=True, match_pseudos=False,
        permutations=[Permutation.Identity, Permutation.Reverse]):
    """Generate a :class:`Forcefield.Plugin` matching atom tuples to a table.

    This is a convenience function to help in the creation of a new user-defined
    forcefield plugin, if the plugin simply matches a particular type of atom
    tuple to a particular param table and writes the matches to a term table in
    the system. Plugins with more complex behavior (e.g. modifying atom mass or
    charge, generating a new parameter table on-the-fly, allowing for multiple
    matches for a single atom tuple, etc.) must be hand-written using the
    :class:`ParameterMatcher` and associated classes.

    'name' should correspond to the name of the plugin as listed in the rules
    file, the name of the parameter file, and the name of the 
    :class:`msys.TermTable` to be added to the system. 'tuple_type' determines 
    what to match; supported values are 'atoms', 'bonds', 'angles', 'dihedrals',
    and 'impropers'. When 'tuple_type' is 'atoms' or 'bonds', 'match_pseudos' 
    determines whether to match virtual/drude particles and bonds from 
    virtuals/drudes to their host atoms. 'permutations' determines which 
    permutations of the atom tuple are allowed during matching, and 'required' 
    determines whether unmatched atom tuples should be ignored or should throw 
    an error.

    'improper' tuples defined by the forcefield are assumed to place the center
    atom first.

    Arguments:
        name -- str

        tuple_type -- str (one of 'atoms', 'bonds', 'angles', 'dihedrals', or 'impropers')

        required -- bool

        match_pseudos -- bool

        permutations -- [f(:class:`Pattern`) -> :class:`Pattern`, ... ]

    Returns: :class:`Forcefield.Plugin`
    """
    def f(tsystem, ff):
        if len(ff.params(name)) == 0:
            raise RuntimeError("'" + name +
                    "' plugin found without '" + name + "' table")
        if tuple_type == 'atoms':
            natoms = 1
            tuples = tsystem.typedAtoms
            if not match_pseudos:
                def lamb(a): return a.atomic_number > 0
                tuples = list(filter(lamb, tuples))
        elif tuple_type == 'bonds':
            natoms = 2
            tuples = tsystem.nonPseudoBonds
            if match_pseudos:
                tuples += tsystem.pseudoBonds
        elif tuple_type == 'angles':
            natoms = 3
            tuples = tsystem.angles
        elif tuple_type == 'dihedrals':
            natoms = 4
            tuples = tsystem.dihedrals
        elif tuple_type == 'impropers':
            natoms = 4
            tuples = tsystem.impropers
        else:
            raise RuntimeError('Unsupported tuple type: '
                    + tuple_type)
        if len(tuples) == 0:
            return
        if tuple_type == 'impropers':
            sys_to_pattern = SystemToPattern.BondToFirst
        else:
            sys_to_pattern = SystemToPattern.Bonded
        matcher = ParameterMatcher.FromFF(ff, name, sys_to_pattern,
                TypeToPattern.Default, permutations)
        table = tsystem.system.addTable(name, natoms,
                Forcefield.ParamTable(name))
        table.category = 'bond'
        for tuple in tuples:
            param, perm = matcher.match(tsystem, tuple)
            if param.id == -1:
                if required:
                    raise RuntimeError('No match found for table '
                            + name + ', pattern '
                            + sys_to_pattern(tsystem, tuple).__str__())
                else:
                    continue
            table.addTerm(tuple, param)
    return Forcefield.Plugin(f)

def PrintParams(params):
    """Print an :class:`msys.ParamTable` to a string.

    Arguments:
        params -- :class:`msys.ParamTable`

    Returns: str

    """
    props = params[0].table.props
    widths = [len(prop) for prop in props]
    table = [props]
    for param in params:
        row = []
        for i, prop in enumerate(props):
            row.append(str(param[prop]))
            widths[i] = max(len(row[-1]), widths[i])
        table.append(row)
    output = ''
    for row in table:
        for i in range(len(props)):
            output += row[i].rjust(widths[i]+1)
        output += '\n'
    return output

def GetBondsAnglesDihedrals(system):
    ''' Return bonds, angles and dihedrals deduced from bond topology
    Returns: dict
    '''
    return _viparr.GetBondsAnglesDihedrals(system._ptr)

def SystemToDot(system, residue_id = msys.BadId):
    ''' Return dot file representation of system
    
    If residue_id is not None, draw atoms not in that residue as external.
    '''
    return _viparr.SystemToDot(system._ptr, residue_id)

def ApplyLigandForcefields(system, ligands, selection='all',
                           rename_atoms=False, rename_residues=False,
                           exhaustive_matching=False,
                           match_bond_stereo=True,
                           match_tet_stereo=False,
                           match_hydrogen=False,
                           verbose=False):
    ''' Apply forcefields from ligands to system

    Args:
        system (msys.System) partially parameterized system
        ligands (list[msys.System]) systems containing one molecule
        selection (str): atom selection for ligands in system
        match_bond_stereo: require matching double-bond stereo in inchi
        match_tet_stereo: require matching tetrahedral stereo in inchi
        match_topology: require matching the hydrogen layer in inchi

    Returns:
        msys.System with new forcefield

    Notes:
        An exception is raised unless each molecule in the selection is
        matched by exactly one ligand.

        Since the parameterized ligand may contain virtual sites not present
        in the input system, no guarantees can be made about preserving the
        atom order, although a best effort is made.

        The match_hydrogen option turns on matching on the hydrogen layer.
        This shouln't be necessary since we already do a Graph match to
        determine a true match between template and ligand, and leaving
        it off by default avoids problems with, e.g. charge imidazole
        having an ambiguous formal charge assignment.
    '''
    def get_inchi(mol, sel='all'):
        lig = mol.clone('(%s) and atomicnumber > 0' % sel)
        inchi = msys.InChI(lig, SNon=False).string

        if not match_hydrogen:
            layers = inchi.split('/')
            inchi = '/'.join(x for x in layers if x[0] not in 'h')

        if not match_bond_stereo:
            layers = inchi.split('/')
            inchi = '/'.join(x for x in layers if x[0] not in 'b')

        if not match_tet_stereo:
            layers = inchi.split('/')
            inchi = '/'.join(x for x in layers if x[0] not in 'tms')

        return inchi

    msys.AssignBondOrderAndFormalCharge(system)

    # construct and validate the selection
    selected = system.select(selection)
    if not selected:
        print("WARNING: Skipping ApplyLigandForcefields -- No atoms in ligand selection '%s'" % selection)
        return system

    fragids = {a.fragid for a in system.select(selection)}
    fragsel = [a for a in system.atoms if a.fragid in fragids]
    if selected != fragsel:
        raise ValueError("Selection did not cover entire molecules: '%s'" % selection)

    # Map ligands by inchi, checking for duplicates
    inchi_to_ligand = dict()
    for i, ligand in enumerate(ligands):
        ligand = ligand.clone()
        msys.AssignBondOrderAndFormalCharge(ligand)
        if {a.fragid for a in ligand.atoms} != {0}:
            raise ValueError("template ligand %s has multiple fragments" % ligand.name)
        inchi = get_inchi(ligand)
        if inchi in inchi_to_ligand:
            raise ValueError("Duplicate inchis '%s' in ligands: %s, %s" % (
                inchi, ligand.name, inchi_to_ligand[inchi].name))
        inchi_to_ligand[inchi] = ligand
        if verbose:
            print("\nLigand %d: %s" % (i, inchi))

    # Assign a ligand to each molecule in the selection
    mappings = list() # [atoms], ligand, inchi, fields
    for fragid in sorted(fragids):
        inchi = get_inchi(system, 'fragid %d' % fragid)
        ligand = inchi_to_ligand.get(inchi)
        if ligand is None:
            raise ValueError("No ligand for molecule with fragid %d inchi %s" % (fragid, inchi))
        atoms = system.select('fragid %d' % fragid)
        fields = dict()
        ct_ids = set(a.residue.chain.ct.id for a in atoms)
        for ct_id in ct_ids:
            ct = system.ct(ct_id)
            for key in ct.keys():
                fields[key] = ct[key]
        mappings.append([atoms, ligand, inchi, fields])
        if verbose:
            print("Matched ligand with fragid %d to ligand %s" % (fragid, ligand.name))

    # for each mapping, construct a graph match between real atoms.
    for mapping in mappings:
        atoms, ligand, inchi, fields = mapping
        mol_atoms = [a for a in atoms if a.atomic_number > 0]
        lig_atoms = [a for a in ligand.atoms if a.atomic_number > 0]
        mol_graph = msys.Graph(mol_atoms)
        lig_graph = msys.Graph(lig_atoms)
        if exhaustive_matching:
            matches = mol_graph.matchAll(lig_graph)
        else:
            match = mol_graph.match(lig_graph)
            matches = [match] if match else []

        # for each match, copy the mol coordinates and make sure the inchi
        # still matches.
        for match in matches:
            for mol_atom, lig_atom in match.items():
                lig_atom.pos = mol_atom.pos
            for a in ligand.atoms:
                if a.atomic_number == 0 and a.nbonds > 0:
                    a.pos = a.bonded_atoms[0].pos
            lig_inchi = get_inchi(ligand)
            if lig_inchi == inchi:
                # inchi still matches, reorder atoms within ligand to match system order
                ordered_atoms=sorted(match.items(), key=lambda x: x[0].id)
                order=[]
                virts=[]
                for mol_atom, lig_atom in ordered_atoms:
                    if(not rename_atoms): lig_atom.name = mol_atom.name
                    # FIXME: This probably isnt the best way to set these
                    if(not rename_residues): 
                        lig_atom.residue.name = mol_atom.residue.name
                    lig_atom.residue.resid     = mol_atom.residue.resid 
                    lig_atom.residue.insertion = mol_atom.residue.insertion
                    lig_atom.residue.chain.name = mol_atom.residue.chain.name
                    lig_atom.residue.chain.segid = mol_atom.residue.chain.segid
                    
                    order.append(lig_atom)
                    # add in virtual atoms at the end
                    virts.extend(sorted([ v for v in lig_atom.bonded_atoms if v.atomic_number == 0]))
                # there shouldnt be duplicated virtuals
                assert(len(virts) == len(set(virts)))
                order.extend(virts)
                # we better have found all the atoms
                assert(len(order) == ligand.natoms)
                # keep the re-ordered ligand
                mapping[1] = ligand.clone(order)
                for k,v in fields.items():
                    mapping[1].ct(0)[k] = v
                break
        else:
            raise RuntimeError("No graph match found for ligand with fragid %d that preserves the inchi %s" % (
                atoms[0].fragid, inchi))

    # now construct a new system, making a best effort to preserve the
    # original atom order.
    mappings.sort(key=lambda x: x[0][0].id)
    ids = { a.id for a in system.atoms }
    newmol = system.clone('none')
    system_fragids = [a.fragid for a in system.atoms]

    for atoms, ligand, inchi, fields in mappings:
        lig_start = atoms[0].id
        initial = { i for i in ids if i < lig_start }
        # preserve bonds!
        initial_fragids = { system_fragids[i] for i in initial }
        initial.update(i for i in ids if system_fragids[i] in initial_fragids)
        if initial:
            newmol.append(system.clone(initial))
            ids.difference_update(initial)
        newmol.append(ligand)
        ids.difference_update(a.id for a in atoms)
    if ids:
        # handle stragglers
        newmol.append(system.clone(ids))

    return newmol
