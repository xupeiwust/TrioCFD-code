"""
This module mainly create an object called :class:`~commons.Tools.Field` which inherit from ``np.ndarray``.
It possesses more features adapted to the use in fluid dynamics :

* the mesh is embbeded in the Field
* the ``Field`` has a name and register the history of the operations leading to its creation
* the ``Field`` possess a method :meth:`~commons.Tools.Field.plot`

"""

# import matplotlib
# import tkinter
# matplotlib.use('TkAgg')
from matplotlib import rc

rc("text", usetex=True)
rc("font", size=14)
rc("font", family="serif")
rc("legend", fontsize=13)
rc("figure", max_open_warning=50)
# rc('figure', figsize=(5, 3.33))
rc("figure", dpi=200)
rc("savefig", dpi=300)
compteur = 0

# from matplotlib.lines import Line2D
# from math import floor, ceil
# from matplotlib import rc
from itertools import product
import numpy as np
import scipy.interpolate as inter

np.set_printoptions(threshold=100)
# from cmath import sin
import matplotlib.pyplot as plt

# import math
import pickle

# from sympy import sympify, latex
from IPython.display import display_latex
import re as re
from itertools import product

# MODIF GAB
try:
    import MEDLoader as ml
except ModuleNotFoundError:
    print("pas passé1")

try:
    import medcoupling as mc
except ModuleNotFoundError:
    print("pas passé2")

try:
    import commons.Bubble_and_Format_Tools as obfrm
except ModuleNotFoundError:
    print("Tools failed to import Bubble_and_Format_Tools")

import time  # Pour la création d'un nouveau maillage irrégulier

# FIN MODIF GAB

font = {"family": "normal", "weight": "bold", "size": 4}
# rc('font', **font)

# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.unicode'] = True
# matplotlib.rcParams["backend"] = "Qt4Agg"
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{gensymb}",r"\usepackage{txfonts}",r"\usepackage{nicefrac}",r"\usepackage{amssymb}",r"\usepackage{sistyle}"]
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{gensymb}",r"\usepackage{txfonts}",r"\usepackage{nicefrac}",r"\usepackage{amssymb}"] #,r"\usepackage{sistyle}"]
# rc('legend', fontsize='medium',numpoints=2)
# rc('text', usetex = True)
# rc('text', usetex = True)
# rc('text', usetex=True)
# rc('font', family='serif')
# import matplotlib.pyplot as pl
spaceDims = 3
dims = range(spaceDims)
HANDLED_FUNCTIONS = {}


def implements(numpy_function):
    """Register an __array_function__ implementation for FieldNew objects."""

    def decorator(func):
        HANDLED_FUNCTIONS[numpy_function] = func
        return func

    return decorator


@implements(np.mean)
def mean(array, *args, axis=-1, **kwargs):
    """

    Args:
        axis (int): the axis along witch to perform the mean
        array (FieldNew): the array to apply the mean operation

    Returns:
        FieldNew: the mean of array (along last axis by default)

    """
    resultat = np.mean(array.view(np.ndarray), axis=axis)
    resultat = array.createField(resultat, name=array.name)
    if resultat.tex is not None:
        resultat.tex = " \\langle " + resultat.tex + " \\rangle "
    if resultat.name is not None:
        resultat.name = " \\langle " + resultat.name + " \\rangle "
    if resultat.label is not None:
        for i, lab in enumerate(resultat.label):
            resultat.label[i] = " \\langle " + lab + " \\rangle "
    return resultat


@implements(np.max)
def max(array, *args, axis=-1, **kwargs):
    """

    Args:
        axis (int): the axis along witch to perform the mean
        array (FieldNew): the array to apply the mean operation

    Returns:
        FieldNew: the max of array (along last axis by default)

    """
    resultat = np.max(array.view(np.ndarray), axis=axis)
    resultat = array.createField(resultat)
    if resultat.tex is not None:
        resultat.tex = "{" + resultat.tex + "}_{max}"
    if resultat.name is not None:
        resultat.name = "{" + resultat.name + "}_{max}"
    return resultat


@implements(np.min)
def min(array, *args, axis=-1, **kwargs):
    """

    Args:
        axis (int): the axis along witch to perform the mean
        array (FieldNew): the array to apply the mean operation

    Returns:
        FieldNew: the min of array (along last axis by default)

    """
    resultat = np.min(array.view(np.ndarray), axis=axis)
    resultat = array.createField(resultat)
    if resultat.tex is not None:
        resultat.tex = "{" + resultat.tex + "}_{min}"
    if resultat.name is not None:
        resultat.name = "{" + resultat.name + "}_{min}"
    return resultat


@implements(np.convolve)
def convolve(field, n_conv, *args, **kwargs):
    """

    Args:
        axis (int): the axis along witch to perform the mean
        array (FieldNew): the array to apply the mean operation

    Returns:
        FieldNew: the min of array (along last axis by default)

    """
    resultat = np.zeros_like(field)
    resultat.name = field.name
    if resultat.tex is not None:
        resultat.tex = "{" + field.tex + "}^{c}"
    for inds in zip(*[range(n) for n in field.shape[:-1]]):
        resultat[inds] = np.convolve(
            field[inds].view(np.ndarray), np.ones((n_conv,)) / n_conv, mode="same"
        )
    return resultat


@implements(np.einsum)
def einsum(op, *args, name=None):
    """
    This method is based on ``np.einsum`` but handles the name and tex arguments

    Args:
        op (str): the ``np.einsum`` operation to be executed
        *args: the list of tensor to give to ``np.einsum``
        name: the name of the result

    Returns:
        Field: a Field with the appropriate value and tex and name

    """
    largs = [arg.view(np.ndarray) for arg in args]
    prod = np.einsum(op, *largs)
    op = op.replace(" ", "")
    op = op.replace("...", "")
    op = op.replace("x", "")
    argsop, resop = op.split("->")
    argsoplist = argsop.split(",")
    prod = args[0].createField(prod)
    prod.label = np.array([])
    prod.tex = ""
    sym = ""
    i = 0
    try:
        for i, (indices, arg) in enumerate(zip(argsoplist, args)):
            # prod.tex += sym + '(' + arg.name + ')_' + indices + ' '
            prod.tex += sym + "{" + arg.name + "}_{" + indices + "} "
            if sym is "":
                sym = " "
    except TypeError:
        print(
            "Warning : The argument number {} has no name, so it can't be added indices".format(
                i
            )
        )
    if name is not None:
        # if resop is not '':
        #     prod.name = '{' + name + '}_{' + resop + '}'
        # else:
        #     prod.name = name
        prod.name = name
    return prod


# TODO: faire en sorte que la méthode gère les vecteurs, tenseurs, et non seulement les champs scalaires
@implements(np.gradient)
def gradient(a, name=None, prints=False, **kwargs):
    """
    This method is based on ``np.gradient`` but handles the name and tex arguments.
    Contrary to the np.gradient methods, it handles unordered structures (such as Fields),
    but in exchange, the computational cost is high.

    Args:
        a (Field): the Field to use gradient on
        name (str): the name of the result
        prints (bool): option de debug pour afficher les résultats intermédiaires

    Returns:
        Field: a Field with the appropriate value and tex and name

    """
    # 1. calcul de la grille ordonnée adaptée en suppossant que les points sont répartis uniformément dans le domaine

    # calcul des deltas entre points pour chaque direction
    grids = []
    for i in range(a._ax.shape[0]):
        xmin = np.min(a._ax[i, :])
        xmax = np.max(a._ax[i, :])
        dx = xmax - xmin
        d = np.min(a._ax[i, :]) - np.mean(a._ax[i, :])
        delta = 6 * (np.std(a._ax[i, :]) ** 2 - d * (d + dx)) / dx - 2 * dx
        grids.append(
            np.linspace(xmin, xmax, int(dx / delta) + 5)
        )  # le +5 est non nécessaire, un +2 aurait été suffisant, mais c'est un moyen d'avoir des grilles un tout petit peu plus fines
    if prints:
        print("grids  ", grids)

    # création de la grille
    mesh = np.meshgrid(*grids, indexing="ij")  # type: tuple[np.ndarray]
    if prints:
        print("mesh ", mesh)

    # 2. Interpolation de a sur la nouvelle grille
    if prints:
        print("ax shape ", np.squeeze(a._ax).shape)
        print("a shape ", np.squeeze(a).shape)
        print("mesh shape ", mesh[0].shape)
    if len(np.squeeze(a).shape) == 1:
        a_interp = inter.griddata(
            np.squeeze(a._ax.T), np.squeeze(a), tuple(mesh), method="cubic"
        )
    else:
        a_interp = np.empty(a.shape[:-1] + mesh[0].shape)
        ranges = [range(a.shape[i]) for i in range(len(a.shape) - 1)]
        for tup in product(*ranges):
            a_interp[tup] = inter.griddata(
                np.squeeze(a._ax.T), a[tup], tuple(mesh), method="cubic"
            )
        # return NotImplementedError
        # for tup in zip(*):
        # TODO: à finir

    if prints:
        print("a_interp ", a_interp)

    # 3. Calcul du gradient
    grads_interp = np.gradient(a_interp, *grids, axis=-1)
    if type(grads_interp) is np.ndarray:
        grads_interp = [grads_interp]
    if prints:
        print("grads_interp ", grads_interp)
        print("grad_interp shape :", grads_interp[0].shape)

    # 4. Interpolation inverse (depuis la grille vers le maillage destructuré)
    # On suppose que tous les grads_interp ont les mm dimensions
    grad = np.empty(a.shape)
    ranges = [
        range(grads_interp[0].shape[i]) for i in range(len(grads_interp[0].shape) - 1)
    ]
    if prints:
        print('ranges : ')
        print(ranges)
    for tup in product(*ranges):
        if prints:
            print('tup : ')
            print(tup)
        grad[tup] = np.squeeze(
            np.array(
                [
                    inter.interpn(tuple(grids), grad_interp[tup], np.squeeze(a._ax.T))
                    for grad_interp in grads_interp
                ]
            )
        )
    # grad = np.squeeze(np.array([inter.interpn(tuple(grids), grad_interp, np.squeeze(a._ax.T)) for grad_interp in grads_interp]))
    if prints:
        print("grad shape", grad.shape)
        print("grad ", grad)

    res = a.createField(grad, name)
    if prints:
        res.print(field=True)
    if a.name is not None:
        res.tex = "\\nabla{" + a.name + "}"
    elif a.tex is not None:
        res.tex = r"\nabla{" + a.tex + "}"
    else:
        res.tex = "grad"
    if res.name is None:
        res.name = res.tex
    if prints:
        res.print(field=True)

    return res


class Field(np.ndarray):
    """
    Cette classe représente un champ nommé porté par un maillage. Elle sub-class numpy ndarray pour tirer directement
    profit des méthodes numpy puissantes. Elle se contente d'ajouter des coordonnées et des labels, et permet la bonne
    gestion de ceux-ci. La représentation du champs vis-à-vis du maillage est non structurée, en effet, pour chaque
    points du maillage, on a le champs qui correspond. Notre array est donc de dimension n_points_maillage x dim_tenseur

    Attributes:
        name (str): the name of the field
        _ax (np.ndarray): the coordinates of the points (the mesh)
        labelRef (str): the label of the coordinates
        label (np.ndarray): the name of each compo
        tex (str): the tex legend corresponding to the Field

    """

    __TRUST_FILES_CACHE = {}

    def __new__(
        cls,
        input_array,
        *args,
        lata_name=None,
        name=None,
        ax=None,
        label_ref=None,
        tex=None,
        label=None,
        **kwargs
    ):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.name = name  # GAB On va dire que c'est le beau nom en mode LaTex
        obj._ax = ax
        obj.labelRef = (
            label_ref  # GAB On va dire que c'est le NomCompo de beau_final.py ...
        )
        obj.tex = tex
        if label is None:
            obj.label = np.array([])
        else:
            obj.label = label
        obj.indices_sel = np.array([None] * len(obj.shape))

        obj.lataName = (
            lata_name  # GAB On va dire que c'est le NomChamp de beau_final.py ...
        )
        obj.listeTemps = []
        obj.mesh = 0
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.name = getattr(obj, "name", None)
        ax = getattr(obj, "_ax", None)
        if ax is None:
            self._ax = None
        else:
            self._ax = ax.copy()
        self.labelRef = getattr(obj, "labelRef", None)
        self.tex = getattr(obj, "tex", None)
        self.label = getattr(obj, "label", np.array([])).copy()
        self.indices_sel = getattr(
            obj, "indices_sel", np.array([None] * len(self.shape))
        ).copy()

    def __getitem__(self, key):
        new = self.createField(super(Field, self).__getitem__(key))
        tot = ""
        deb = r"\left["
        fin = r"\right]"
        if type(key) is not tuple:
            key = (key,)
        for i, sl in enumerate(key):
            if type(sl) is slice:
                if sl.start is None:
                    start = ""
                else:
                    start = str(sl.start)
                if sl.stop is None:
                    stop = ""
                else:
                    stop = str(sl.stop)
                if (sl.step is None) or (sl.step == 1):
                    step = ""
                else:
                    step = str(sl.step)
                if start == stop == step == "":
                    ind = ""
                else:
                    ind = start + ":" + stop + ":" + step
            elif sl is None:
                ind = ""
            else:
                ind = str(sl)
            tot += ind
            if i != (len(key) - 1):
                tot += r",\ "
        if tot.replace(r",\ ", "") != "":
            tot = deb + tot + fin

        new.indices_sel = np.array(key)
        if self.name is not None:
            new.name = self.name
            new.name = new.name + tot
        if self.tex is not None:
            new.tex += tot
        return new

    def __array_ufunc__(self, ufunc, method, *args, sym=" ? ", **kwargs):
        """
        Overrides all of the ufunctions (add, mul, sub, ...)
        Create a new array by applying the original ufunc called, and set it tex variable in accordance to the operation

        >>> a = Field([0,0], name='a')
        >>> b = Field([0,1], name='b')
        >>> c = a + b
        >>> c.tex
        '(a + b)'

        Args:
            ufunc:
            method:
            *args:
            sym:
            **kwargs:

        Returns:

        """
        args_array = []
        for arg in args:
            # if isinstance(arg, Field) or isinstance(arg, np.ndarray) or isinstance(arg, float) \
            #         or isinstance(arg, int):
            args_array.append(np.array(arg))
        # print(args_array)
        # print(ufunc)
        # print(args)
        res = super(Field, self).__array_ufunc__(ufunc, method, *args_array, **kwargs)
        # print(res)

        if res is NotImplemented:
            return NotImplemented

        res = Field(res, ax=self._ax, label_ref=self.labelRef)
        # res = self.createField(res, name=self.name)

        # Checking the ax compatibility
        if res._ax is None:
            for arg in args:
                if isinstance(arg, Field):
                    if arg._ax is not None:
                        res._ax = arg._ax
                        break

        if res._ax is not None:
            for arg in args:
                if isinstance(arg, Field):
                    if arg._ax is not None:
                        if not np.all(res._ax == arg._ax):
                            if arg.name is not None:
                                name0 = arg.name
                            else:
                                name0 = str(arg)
                            if self.name is not None:
                                name1 = self.name
                            else:
                                name1 = "self"
                            diff = np.max(np.abs(res._ax - arg._ax))
                            tol = 1e-12
                            if diff > tol:
                                print(
                                    "Warning : {} and {} doesn't have the same ax !!!".format(
                                        name0, name1
                                    )
                                )
                                print(
                                    "Difference=%g is larger than tol=%g" % (diff, tol)
                                )

        f = ufunc.__name__
        # print(f)
        # li = []
        # for arg in args:
        #     try:
        #         if arg.name is not None:
        #             li.append(arg.name)
        #         else:
        #             li.append(arg.tex)
        #     except:
        #         li.append(None)
        # print(li)
        # print('self', self)
        # print('arg0', args[0])

        # préfixe de l'opération

        res.tex = ""
        start_with_parenth = False
        if len(args) == 1:
            start_with_parenth = True
            if f == "sqrt":
                res.tex = "\\sqrt{"
            elif f == "sin":
                res.tex = "\\sin{"
            elif f == "cos":
                res.tex = "\\cos{"
            elif f == "tan":
                res.tex = "\\tan{"
            elif f == "square":
                args += (2,)
                start_with_parenth = False
            else:
                # res.tex = '('
                start_with_parenth = False
        if len(args) == 2:
            if f == "add":
                sym = " + "
            elif f == "subtract":
                sym = " - "
            elif f == "multiply":
                sym = " * "
            elif f == "power" or f == "square":
                sym = "}^{"
            elif f == "divide" or f == "true_divide":
                sym = " / "
            else:
                sym = " ? "
                # print(sym)
        if len(args) > 2:
            raise Exception("Too much arguments in this operation, not implemented yet")

        arg = args[0]
        try:
            if arg.name is not None:
                res.tex += arg.name
            else:
                res.tex += arg.tex
        except:
            res.tex += str(arg)

        # if sym == ' ? ':
        #     print(res.tex)

        # deuxième membre

        # au milieu de l'opération avec le deuxième argument (s'il y en a un)
        for arg in args[1:]:
            # print(arg)
            try:
                if arg.labelRef is not None:
                    res.labelRef = arg.labelRef
                if arg.name is not None:
                    # res.tex += sym + arg.name
                    res.tex = addpar(res.tex, arg.name, sym)
                else:
                    # res.tex += sym + arg.tex
                    res.tex = addpar(res.tex, arg.tex, sym)
            except:
                # res.tex += sym + str(arg)
                res.tex = addpar(res.tex, str(arg), sym)

        # pour finir
        if start_with_parenth:
            res.tex += "}"

        return res

    def __array_function__(self, func, types, *args, **kwargs):
        """
        This mehtod allows to override numpy functions with functions decorated with :py:func:`implements`

        Args:
            func:
            types:
            args:
            kwargs:

        Returns:

        """
        if func not in HANDLED_FUNCTIONS:
            return super(Field, self).__array_function__(func, types, *args, **kwargs)
        # Note: this allows subclasses that don't override
        # __array_function__ to handle MyArray objects
        if not all(issubclass(t, type(self)) for t in types):
            return NotImplemented
        return HANDLED_FUNCTIONS[func](*args, **kwargs)

    # TODO: gérer l'affichage de maillage 2D et plus
    def get_tex(self, p=False, val=False, name=True, tex=True, without_dollar=False):
        """
        Returns:
            str: the latex text corresponding to the operations strored in self.tex

        Examples:
            >>> a = Field([-1, 1], name='a')
            >>> b = Field([1, -1], name='b')
            >>> c = a + b
            >>> c = a * np.sqrt(c)
            >>> c.get_tex()
            'a \\sqrt{a + b}'

        """
        members = []
        if name and self.name is not None:
            try:
                # members.append(latex(sympify(self.name, evaluate=False)))
                members.append(self.name)
            except:
                members.append("")
        elif name and (not tex) and (self.tex is not None):
            try:
                # members.append(latex(sympify(self.tex, evaluate=False)))
                members.append(self.tex)
                tex = False
            except:
                members.append("")
        if tex and self.tex is not None:
            try:
                # members.append(latex(sympify(self.tex, evaluate=False)))
                members.append(self.tex)
            except:
                members.append("")
        if val:
            if (self.size == 1) and not (self.dtype == "<U9"):
                members.append("%.3g" % self)
            elif not (self.dtype == "<U9"):
                try:
                    # members.append(latex(sympify(str(self), evaluate=False)))
                    members.append(str(self))
                except:
                    members.append("")
            else:
                members.append(r"\textrm{%s}" % self)

        if len(members) == 0:
            return ""
        la = ""
        for mem in members[:-1]:
            la += mem + " = "
        la += members[-1]
        # n = 1
        # la = re.sub(r'\\operatorname{mean}{\\left\((.*?)\\right\)}', r'\\langle \1 \\rangle', la)
        # la = re.sub(r'operatorname{(.*?)}{\\left\((.*?)\\right\)}', r'\1{\2}', la)
        # la = re.sub(r'operatorname{(.*?)}\^{(.*?)}{\\left\((.*?)\\right\)}', r'\1\^{\2}{\3}', la)
        # la = re.sub(r'operatorname{(.*?)}_{(.*?)}{\\left\((.*?)\\right\)}', r'\1\^{\2}{\3}', la)
        # la = re.sub(r"_{p}", "'", la)
        if not without_dollar:
            la = "$" + la + "$"
        if p:
            if without_dollar:
                la = "$" + la + "$"
            display_latex(la, raw=True)
        return la

    def to_dict(self):
        if self._ax is not None:
            ax = self._ax.tolist()
        else:
            ax = None
        if isinstance(self.label, list):
            label = self.label
        else:
            label = self.label.tolist()
        res = {
            "value": self.tolist(),
            "name": self.name,
            "tex": self.tex,
            "ax": ax,
            "label_ref": self.labelRef,
            "lataName": self.lataName,
            "label": label,
        }
        return res

    @classmethod
    def from_dict(cls, dic):
        ax = dic["ax"]
        if ax is not None:
            ax = np.array(ax)
        return cls(
            np.array(dic["value"]),
            ax=ax,
            lata_name=dic["lataName"],
            name=dic["name"],
            tex=dic["tex"],
            label_ref=dic["label_ref"],
            label=np.array(dic["label"]),
        )

    def print(self, field=False, **kwargs):
        """
        Present the Field in a nice print with all attributes if field is True.
        Otherwise only prints the name and tex of the Field.

        Examples:

            .. jupyter-execute::

               a = Field(np.arange(9*5).reshape((3,3,5)), name='A',
                          ax=np.arange(5)/5)
               a.print()

        """
        display_latex(self.get_tex(**kwargs), raw=True)
        if field:
            print("Tenseur de taille : ", self.shape)
            value = str(self)
            for i in range(len(self.shape), 0):
                value.replace("[" * i, " " * (len(self.shape) - i) + "[" * i)
            print("    " + value.replace("\n", "\n    "))
            for key, value in sorted(self.__dict__.items()):
                if key == "_ax":
                    value = str(value).replace("\n", "\n           ")
                print("    " + key + " : ", value)
            if len(self.shape) > 0:
                me = mean(self)
                value_me = str(me)
                for i in range(len(self.shape), 0):
                    value_me.replace("[" * i, " " * (len(me.shape) - i) + "[" * i)
                print("    mean : " + value_me.replace("\n", "\n           "))

    def createField(self, array=None, name=None):
        """
        creation of a Field from an other Field

        Args:
            name (str): a name for the new field

        """
        if array is not None:
            other = Field(
                array, name=name, ax=self._ax, label_ref=self.labelRef, tex=self.tex
            )
        else:
            other = Field(
                self[:].copy(),
                name=name,
                ax=self._ax,
                label_ref=self.labelRef,
                tex=self.tex,
            )
        return other

    def grad(self):
        return np.gradient(self)

    @property
    def _npa(self):
        """
        For back compatibility with Field_old

        Returns:
            np.ndarray: self[...]

        """
        return self.view(np.ndarray)

    @_npa.setter
    def _npa(self, array):
        self[...] = self.createField(array)

    def settex(self, tex):
        self.tex = tex

    def setname(self, name):
        self.name = name

    def settexname(self, txt):
        self.name = txt
        self.tex = r"$%s$" % txt

    # DONE: il faudra supprimer les parenthèses dans l'usage de TypeField
    @property
    def TypeField(self):
        """
        Return the type of the object
        """
        # DONE: remplacer isinstance par condition sur len self
        shape = self.shape[:-1]
        if shape == (1,):
            a = "scalar"
        elif len(shape) == 1:
            a = "vector"
        else:
            a = "tensor{}".format(len(shape))
        return a

    # def grad(self, bc=0):
    #     """
    #     Return the gradient of self

    #     Args:
    #          bc (int): the boundary conditions type

    #     Returns:
    #         the gradient of the object

    #     """
    #     # TODO: supprimer la distinction des cas ou utiliser dans les différentes classes

    #     other = self.derivate(bc)
    #     a = self.TypeField

    #     if a == 'scalar':
    #         aa = len(self[0, :])
    #         grad = other.createField(np.zeros((3, aa)))
    #         grad[2, :] = other[0, :]

    #     elif a == 'vector':
    #         # TODO: pourquoi on recalcule derivate ?
    #         # si c'est inutile, on peut tout faire sans distinctions
    #         other = self.derivate()
    #         c = len(other._ax[0, :])
    #         grad = other.createField(np.zeros((3, 3, c)))
    #         grad[:, 2, :] = other

    #     else:
    #         c = other._ax.shape[1]
    #         shape = [3] * len(self.shape)
    #         shape += [c]
    #         grad = other.createField(np.zeros(shape))
    #         grad[..., 2, :] = other

    #     if self.name is not None:
    #         grad.name = 'grad(' + self.name + ')'
    #     if self.tex is not None:
    #         grad.tex = 'grad(' + self.tex + ')'

    #     return grad

    def integ_const(self, inter=None):
        """
        Return the integrated value

        Args:
             inter (list of int): interval

        """

        # TODO: changer ça pour gagner en vitesse d'execution
        dx = self._ax[0, 2] - self._ax[0, 1]
        som = 0.0

        for i in range(len(self[:, 0])):
            for j in range(len(self[i, :])):
                if inter[0] < j < inter[1]:
                    som = som + self[i, j] * dx

        # TODO: gérer l'attribut tex
        return som

    def integ(self):
        """
        Return the integrated field
        """
        inte = np.zeros_like(self)
        if self.name is not None:
            inte.name = "Integral(" + self.name + "," + self.labelRef + ")"
        if self.tex is not None:
            inte.tex = "Integral(" + self.tex + "," + self.labelRef + ")"

        # TODO: boucle
        for i in range(0, len(self[:, 0])):
            somme = 0
            for j in range(0, len(self[i, :]) - 1):
                somme = somme + self[i, j] * (self._ax[0, j + 1] - self._ax[0, j])
                inte[i, j] = somme

        inte[:, -1] = 2 * inte[:, -2] - inte[:, -3]
        return inte

    # TODO: remplacer par convolve
    def MoyenneGlissante(self, db=0.4, maxcurv=0.0, wall=0.1):
        """
        Return the integrated field

         Args:
             db:
             maxcurv (float):
             wall:

         Returns:
             Field: a field of the same dimension as self

        """
        inte = self.copy()
        inte.name = "/[" + self.get_tex(tex=False) + "]d" + self.labelRef
        # inte = self.createField('/[' + self.name + ']d' + self.labelRef)
        cc = 0.0
        for i in range(self.shape[-1]):
            iter = 0
            inte[..., i] = 0.0
            for j in range(self.shape[-1] - 1):
                if (
                    wall < self._ax[0, j] < self._ax[0, i]
                    and self._ax[0, j] < wall + db
                ):
                    cc = cc + 1.0
                    inte[..., i] = inte[..., i] + self[..., j]
                    iter = iter + 1.0

            # if self._ax[0, i] <= wall * 1.3:
            #     inte[..., i] = maxcurv

            if iter != 0:
                inte[..., i] = inte[..., i] / iter

        inte[..., -1] = inte[..., -2]
        return inte

    # TODO: attention la méthode est déja déclarée pour les ndarray

    def mean_compo(self, compo):
        """
        Return the mean value of self

        Args:
             compo: a component

        """
        return mean(self[compo, :])

    def derivate(self, bc_type=0):
        """
        Return the derivative of self

        Args:
            bc_type:

        """
        # ne sert que dans grad pour le moment...
        a = len(self[0, :])
        b = len(self[:, 0])
        other = self.createField(np.zeros((b, a)))
        if self.name is not None:
            self.name = "d(" + self.name + ")/d" + self.labelRef
        if self.tex is not None:
            self.tex = "d(" + self.tex + ")/d" + self.labelRef
        i = 0
        # TODO: ENLEVER CETTE BOUCLE WHILE
        # pour ça il faudrait que je fusionne cell_to_cell_gradient_scalar et cette fonction, ou que je modifie
        # cell_to_cell_... . Je pense qu'il suffit d'indiquer une ellipse dans chaque instance de Field ou autre
        while i <= len(self[:, 0]) - 1:
            other[i, :] = cell_to_cell_gradient_scalar(
                self[i, :], self._ax[0, :], bc_type
            )
            i += 1

        return other

    # TODO: vérifier quelle est la sortie de cette fonction à cause de __array_ufunc__
    def sqrt_pos(self):
        """
        Do not returns a true square root

        Args:
            a Field

        Returns:
            other: the square root of elem, if elem > 0 and -sqrt(-elem) if elem<0 for elem in self

        """
        other = self.createField("sqrt")
        other = np.zeros_like(self)
        pos = self >= 0.0
        neg = self < 0.0
        other[pos] = np.sqrt(self[pos])
        other[neg] = -np.sqrt(-self[neg])
        return other

    # TODO: regarder ce qu'on peut faire au niveau de tex
    def symetriser(self, symet=None):
        """
        Return the symetric of Field
        Attention, ne fait pas l'opération en place

        Args:
            symet (np.ndarray):  un array avec des compo 1 for symetry. -1 for antisymetry. Default is False
        """
        #### utilise

        sym = self.pre_reshape().copy()

        if symet is not None:
            symm = symet > 0
            # symm = symm[..., np.newaxis]
            # symm = np.repeat(symm, self.shape[-1], axis=-1)
            sym[symm, :] = (sym[symm, :] + sym[symm, ::-1]) / 2.0
            antisym = symet < 0
            # antisym = antisym[..., np.newaxis]
            # antisym = np.repeat(antisym, self.shape[-1], axis=-1)
            sym[antisym, :] = (sym[antisym, :] - sym[antisym, ::-1]) / 2.0

        sym = self.post_reshape(sym)
        if self.name is not None:
            sym.name = "{" + self.name + "}_{s}"
        for i in range(len(sym.label)):
            sym.label[i] = "{" + sym.label[i] + "}_{s}"

        return sym

    def pre_reshape(self, arr=None):
        if arr is None:
            arr = self
        s = arr.shape
        if len(s) == 3:
            return arr.reshape((9, -1))
        elif len(s) == 4:
            return arr.reshape((27, -1))
        else:
            return arr

    def post_reshape(self, arr=None):
        if arr is None:
            arr = self
        s = arr.shape
        if s[0] == 9:
            return arr.reshape((3, 3, -1))
        if s[0] == 27:
            return arr.reshape((3, 3, 3, -1))
        else:
            return arr

    def contract(self, ref, indice):
        """
        This method is meant to replace the methods contract product and tensorProduct

        Args:
            ref (int or str): either 0, 'i', 'ij' or 'ijk', is the result wanted after contraction
            indices (str): the operation to compute. Ex: 'ikkj' or 'ii'

        Returns:
            Field: the result of the operation

        """
        if ref == 0:
            ref = ""
        elif len(ref) > 3:
            raise Exception(
                "ref= " + str(ref) + " is not an acceptable entrance, len is too long"
            )
        op = indice + "..." + "->" + ref + "..."
        product = einsum(op, self, name=self.name)
        return product

    def product(self, other, sym="*"):
        """
        This method is for backward compatibility. The prefered way is now to use the operator *
        """
        if sym is "*":
            res = self * other
        elif sym is "/":
            res = self / other
        else:
            raise Exception(
                "The symbol {} is not recognized, symbols recognized are * and /".format(
                    sym
                )
            )

        # verification qu'il y a les mêmes résultats
        if not np.all(res == self.product_old(other, sym)):
            raise Exception("prduct and product_old does not have the same behavior")
        return res

    def MatProduct(self, other):
        """a supprimer a terme par TensorProduct. cette méthode est conservée pour rétro compatibilité"""
        b = other.TypeField
        if b == "tensor2":
            pro = self.product(other)
            aa = len(self[0, :])
            ab = len(self[:, 0])
            res = np.zeros((ab, aa))
            res = self.createField(res, name="not implemented")
            res[0, :] = pro[0, 0, :] + pro[0, 1, :] + pro[0, 2, :]
            res[1, :] = pro[1, 0, :] + pro[1, 1, :] + pro[1, 2, :]
            res[2, :] = pro[2, 0, :] + pro[2, 1, :] + pro[2, 2, :]

        elif b == "vector":
            aa = len(self[0, :])
            ab = len(self[:, 0])
            res = np.zeros((ab, ab, aa))
            res = self.createField(res, name="not implemented")
            for i in dims:
                for j in dims:
                    res[i, j, :] = self[i, :] * other[j, :]
        else:
            raise Exception("NotImplemented")

        return res

    def tensorProduct(self, other, name=None):
        """
        Cette méthode est conservée pour rétro compatibilité, lui préférer l'utilisation de np.einsum(op, a, b)
        """
        if self.TypeField == "vector":
            if other.TypeField == "vector":
                res = einsum("ix,jx->ijx", self, other, name=name)
            elif other.TypeField == "tensor2":
                res = einsum("ix,jkx->ijkx", self, other, name=name)
            else:
                print(
                    "Types in tensor product are: ",
                    self.TypeField,
                    " and ",
                    other.TypeField,
                )
                raise NotImplementedError
        elif self.TypeField == "tensor2":
            if other.TypeField == "vector":
                res = einsum("ijx, kx->ijkx", self, other, name=name)
            elif other.TypeField == "tensor2":
                res = einsum("ijx, klx->ijklx", self, other, name=name)
            else:
                print(
                    "Types in tensor product are: ",
                    self.TypeField,
                    " and ",
                    other.TypeField,
                )
                raise NotImplementedError
        elif self.TypeField == "tensor3":
            if other.TypeField == "vector":
                res = einsum("ijkx, lx->ijklx", self, other, name=name)
            elif other.TypeField == "tensor2":
                res = einsum("ijkx, lmx->ijklmx", self, other, name=name)
            else:
                print(
                    "Types in tensor product are: ",
                    self.TypeField,
                    " and ",
                    other.TypeField,
                )
                raise NotImplementedError
        else:
            print(
                "Types in tensor product are: ",
                self.TypeField,
                " and ",
                other.TypeField,
            )
            raise NotImplementedError

        if name is None:
            res.name = "(" + self.name + ") x (" + other.name + ")"
            res.tex = "(" + self.tex + ") \otimes (" + other.tex + ")"
        return res

    def transpo(self, *args, **kwargs):
        """
        Transpo using einsum
        """
        if len(self.shape) == 2:
            op = "iz -> iz"
        elif len(self.shape) == 3:
            op = "ijz -> jiz"
        else:
            raise NotImplementedError
        # res = np.einsum(op, self, '(' + self.name + ')**T')
        res = einsum(op, self)
        res.name = "(" + self.name + ")**T"
        res.tex = "(" + self.tex + ")^T"
        return res

    def fromScatoVecenz(self):
        """
        Pour le moment méthode créée pour transformer un scalaire issu de la dérivée en z d'un scalaire (type temperature)
        en vecteur 3 avec les 2 premières composantes nulles (gradient_x = 0, gradient_y = 0)
        """
        data = np.zeros((3, self.shape[-1]))
        data[2, :] = self[...]
        return self.createField(data, name=self.name)

    def fromVectoMat33enz(self):
        """
        Pour le moment méthode créée pour transformer un vecteur issu de la dérivée en z d'un vecteur (type vitesse)
        en matrice 3x3 avec les 2 premières composantes nulles (jacobienne)
        """
        data = np.zeros((3, 3, self.shape[-1]))
        data[:, 2, :] = self[...]
        return self.createField(data, name=self.name)

    def fromMat33toMat333enz(self):
        """
        Pour le moment méthode créée pour transformer un mat33 issu de la dérivée en z d'un mat33 (type Rij)
        en matrice 3x3x3 avec les 2 premières composantes (derivee selon x et y) nulles (jacobienne)
        """
        data = np.zeros((3, 3, 3, self.shape[-1]))
        data[:, :, 2, :] = self[...]
        return self.createField(data, name=self.name)

    def fromMat333toMat3333enz(self):
        """
        Pour le moment méthode créée pour transformer un mat333 issu de la dérivée en z d'un mat333 (type ukRij)
        en matrice 3x3x3x3 avec les 2 premières composantes (derivee selon x et y) nulles (jacobienne)
        """
        data = np.zeros((3, 3, 3, 3, self.shape[-1]))
        data[:, :, :, 2, :] = self[...]
        return self.createField(data, name=self.name)

    def ChangerMaillage(self, other):
        """
        Cette methode doit être réécrite
        """
        tmp = other * 0.0
        for i in range(len(tmp._ax[0, :])):
            for j in range(len(self._ax[0, :])):
                if self._ax[0, j] > tmp._ax[0, i]:
                    avant = tmp._ax[0, i] - self._ax[0, j - 1]
                    apres = self._ax[0, j] - tmp._ax[0, i]
                    av = j - 1
                    ap = j
                    break
        print("verif maille", avant + apres, av, ap)
        tmp[:, i] = (avant * self[:, ap] + apres * self[:, av]) / (avant + apres)
        return tmp

    @classmethod
    def LoadFromFile(cls, fileObj, col_ref, cols, name, labelRef, tex, bc=0):
        #### utilise
        yoyo = []
        x, tmp = fileObj.getValues(col_ref[0])
        # print(col_ref[0])
        # print(x)
        # print(tmp)
        yoyo.append(tmp)

        if len(cols) == 1:
            l = []
            # ret = ScalarField(name, np.vstack(yoyo), labelRef, tex)
            x, tmp = fileObj.getValues(cols[0])
            l.append(tmp)
            if bc == 2:
                # TODO: impossible que ça marche, l est toujours de taille 1
                l[0] = (l[1] + l[-2]) * 0.5
                l[-1] = (l[1] + l[-2]) * 0.5
            # print('array : ', np.vstack(l))
            # print('ax : ', np.vstack(yoyo))
            # l'utilisation de l permet de passer d'un array 1d à un array 2d avec 1 seul élément dans la première dir
            ret = cls.__new__(
                cls,
                np.vstack(l),
                name=name,
                ax=np.vstack(yoyo),
                label_ref=labelRef,
                tex=tex,
            )
            # print('arr._ax : ', ret._ax)
            ret.label = cols
            # ret._npa_x = x

        elif len(cols) == spaceDims:
            # ret = VectorField(name, np.vstack(yoyo), labelRef, tex)
            for i in dims:
                x, tmp = fileObj.getValues(cols[i])
                if i == 0:
                    nz = tmp.shape[0]
                    l = np.empty((spaceDims, nz), dtype=np.float64)
                    l[i, :] = tmp
                    if bc == 2:
                        l[i, 0] = (l[i, 1] + l[i, -2]) * 0.5
                        l[i, -1] = (l[i, 1] + l[i, -2]) * 0.5
                else:
                    if tmp.shape[0] != nz:
                        raise Exception
                    l[i, :] = tmp
                    if bc == 2:
                        l[i, 0] = (l[i, 1] + l[i, -2]) * 0.5
                        l[i, -1] = (l[i, 1] + l[i, -2]) * 0.5
                    pass
            ret = cls.__new__(
                cls,
                np.vstack(l),
                ax=np.vstack(yoyo),
                name=name,
                label_ref=labelRef,
                tex=tex,
            )
            # ret._npa_x = x
            ret.label = cols

        elif len(cols) == spaceDims * spaceDims:
            for i in dims:
                for j in dims:
                    x, tmp = fileObj.getValues(cols[i + spaceDims * j])
                    if i == 0 and j == 0:
                        nz = tmp.shape[0]
                        l = np.empty((spaceDims, spaceDims, nz), dtype=np.float64)
                        l[i, j, :] = tmp
                        if bc == 2:
                            l[i, j, 0] = (l[i, j, 1] + l[i, j, -2]) * 0.5
                            l[i, j, -1] = (l[i, j, 1] + l[i, j, -2]) * 0.5
                    else:
                        if tmp.shape[0] != nz:
                            raise Exception
                        l[i, j, :] = tmp
                        if bc == 2:
                            l[i, j, 0] = (l[i, j, 1] + l[i, j, -2]) * 0.5
                            l[i, j, -1] = (l[i, j, 1] + l[i, j, -2]) * 0.5

            ret = cls.__new__(
                cls, l, name=name, label_ref=labelRef, tex=tex, ax=np.vstack(yoyo)
            )
            # ret._npa_x = x
            ret.label = cols

        elif len(cols) == spaceDims * spaceDims * spaceDims:
            # ret = Tensor3Field(name, np.vstack(yoyo), labelRef, tex)
            for i in dims:
                for j in dims:
                    for k in dims:
                        x, tmp = fileObj.getValues(
                            cols[i + spaceDims * j + spaceDims * spaceDims * k]
                        )
                        if i == 0 and j == 0 and k == 0:
                            nz = tmp.shape[0]
                            l = np.empty(
                                (spaceDims, spaceDims, spaceDims, nz), dtype=np.float64
                            )
                            l[i, j, k, :] = tmp
                            if bc == 2:
                                l[i, j, k, 0] = (l[i, j, k, 1] + l[i, j, k, -2]) * 0.5
                                l[i, j, k, -1] = (l[i, j, k, 1] + l[i, j, k, -2]) * 0.5
                        else:
                            if tmp.shape[0] != nz:
                                raise Exception
                            l[i, j, k, :] = tmp
                            if bc == 2:
                                l[i, j, k, 0] = (l[i, j, k, 1] + l[i, j, k, -2]) * 0.5
                                l[i, j, k, -1] = (l[i, j, k, 1] + l[i, j, k, -2]) * 0.5

            ret = cls.__new__(
                cls, l, name=name, ax=np.vstack(yoyo), label_ref=labelRef, tex=tex
            )
            # ret._npa_x = x
            ret.label = cols

        else:
            raise Exception("Number of columns %d given not recognized" % (len(cols)))

        return ret

    ##### MODIF GAB
    # --> A tester une fois qu'on pourra charger des objets 2d et plus
    # ~~~~
    # fmed = '/data/calculsfme/gr262753/'+dossier+'/'+fichier
    # NomRef, NomCompo,NomBeau = "VELOCITY_ELEM_DOM_dual", ["Ux","Uy","Uz"], r'$\bm{u_{SIMU}}$'
    # SIMU_Reb400 = dtool.Simu(jdd);
    # SIMU_Reb400.it0, SIMU_Reb400.itf = 147, 197
    # Reb400_U = SIMU_Reb400.get_field_from_med(fmed,NomRef,NomCompo,NomBeau)
    # ~~~~

    # def fromVectoMat33enz(self):
    #     """
    #     Pour le moment méthode créée pour transformer un vecteur issu de la dérivée en z d'un vecteur (type vitesse)
    #     en matrice 3x3 avec les 2 première composantes nulles (jacobienne)
    #     """
    #     data = np.zeros((3,3,self.shape[-1]))
    #     data[:,2,:] = self[...]
    #     return self.createField(data, name=self.name)

    def LoadFromMED(self, MEDFileName, it0, itf, nb, dual_Mesh=None, give_Mesh=False):
        """
        From a .med file, extracts a field and its timeline

        Args:
             MEDFileName: Nom+(chemin en dur) du fichier .med d'ou extraire le champ
             NomChamp: nom du champ tel qu'il est ecrit dans le .lata
             NomCompo: nom des composantes du champ tels qu'elles apparaîtrons dans les plots et
                       dans Paraview
             NomBeau: nom du champ tels qu'elles apparaîtrons dans les plots et
                       dans Paraview
             labelRef (int): Comme on prend du 3D on s'en sert pas
             tex: ???
             it0, itf: Intervalle des itérations du med sur lesquelles travailler
             jddName: Nom+(chemin en dur) du jeu de données : le .data
             dual_Mesh: Dans le cas ou on recupere un champ sur DOM, il faut aussi donner le maillage
                        du dual pour pouvoir faire l'INTERPOLATION : dans ce cas, dual_Mesh
                        correspond au maillage du dual.
                        Dans le cas ou on récupère un champ sur DOM_dual, laisser dual_Mesh à
                        None.
             give_Mesh: Mettre à True si l'on veut récupérer le maillage
                        Mettre à False sinon (économise de la mémoire)

        Returns:
             NPChamp: [Nx, Ny, Nz, Nc, Nb, Nt], le champ en NumPy (ndarray)
             ListeTemps: [Nt], liste des instants sauvegardés par le code

        """

        # self.name     --> NomChamp, MEDFieldName
        # self.labelRef --> NomCompo
        # self.
        print("in Tools.LoadFromMED : self.name")
        print(self.name)

        ListeMEDFileField, MEDFileMesh, TOF = obfrm.MEDFileToMEDField_SomeTS(
            MEDFileName, self.name, it0, itf, -1
        )
        it = [(i, -1) for i in range(it0, itf)]

        # Maillage
        # if give_Mesh:
        # La c'est qu'on veut récupérer un maillage
        if "dual" in self.name:
            # on ne récupère que le DOM_dual
            dual_Mesh = MEDFileMesh.getMeshAtLevel(0)
            dual_Mesh, _ = obfrm.fromUtoC(
                dual_Mesh
            )  # cet appel prend 8sec pour un maillage 78*32*32

            Mesh_DOM = False
        else:
            # On récupère le DOM qu'on interpolera ensuite
            # sur DOM_dual
            Mesh_DOM = MEDFileMesh.getMeshAtLevel(0)
            Mesh_DOM, _ = obfrm.fromUtoC(
                Mesh_DOM
            )  # cet appel prend 8sec pour un maillage 78*32*32
        del _

        ### En vrai ça, ça doit etre dans le self quand il lit le jdd je pense ...
        xn = dual_Mesh.getCoordsAt(0)
        Nx = len(xn) - 1
        yn = dual_Mesh.getCoordsAt(1)
        Ny = len(yn) - 1
        zn = dual_Mesh.getCoordsAt(2)
        Nz = len(zn) - 1
        ### a voir si on l'a pas direct du jdd plutot ...

        Super_Numpy = np.zeros((Nx, Ny, Nz, nb, len(it)))
        listeTemps = []
        for i, o in it:
            di = i - it0
            o = -1
            t = ml.GetTimeAttachedOnFieldIteration(MEDFileName, self.name, i, o)

            # 1. Récupère le champ de valeurs
            Numpy = ListeMEDFileField[di].getFieldAtLevel(ml.ON_CELLS, 0)
            # S'il se trouve que le champ n'est pas structure, il faudra appliquer en-tete
            if Mesh_DOM == False:
                _, _, _, Numpy = obfrm.MEDFieldToNumpy(Numpy, dual_Mesh)
            else:
                _, _, _, Numpy = obfrm.MEDFieldToNumpy(Numpy, Mesh_DOM)

            if not ("dual" in self.name):
                # pas dual <--> champ aux centres de DOM, à INTERPOLER sur les centres
                # du dual
                # Récupère le champ d'indicatrices
                Numpy.setNature(mc.ExtensiveMaximum)
                rem_Numpy = mc.MEDCouplingRemapper()
                rem_Numpy.prepare(Mesh_DOM, dual_Mesh, "P0P0")
                rem_Numpy = rem_Numpy.transferField(Numpy, 1e300)
                del rem_Numpy
                _, _, _, Numpy = obfrm.MEDFieldToNumpy(rem_Numpy, dual_Mesh)

            Numpy = obfrm.correction3D_I(Numpy)

            if di == 0:
                nc = Numpy.shape[-1]
                Super_Numpy = np.zeros((Nx, Ny, Nz, nc, len(it)))

            Super_Numpy[..., di] = Numpy
            listeTemps.append(t)

        Super_Numpy.transpose(0, 1, 2, 4, 3)
        del (Mesh_DOM, it)
        if give_Mesh:
            self.mesh = dual_Mesh
        del dual_Mesh
        self.NPChamp = Super_Numpy
        del Super_Numpy
        return listeTemps

    # A tester ...
    # ~~~~~
    # post_tt_file = ''
    # BxMesh = SIMU_Reb400.getBxMesh()
    # Reb400_U.SaveToMED(post_tt_file, BxMesh, 'instantanne')
    # ~~~~~
    def SaveToMED(self, fmed_out, BxMesh, temps="moyenne"):
        """
         Writes the field in a .med file.
        Args :
           fmed_out : name of the output file
           BxMesh   : Mesh on which the field lies
           temps    : whether it is a one-frame field or a multiple-frame field.
        """

        MED_FIELD = ml.MEDCouplingFieldDouble.New(0, ml.ONE_TIME)
        MED_FIELD.setName(self.FieldName)
        MED_FIELD.setMesh(BxMesh)

        if temps == "moyenne":
            # On n'écrit qu'un seul pas de temps du coup
            MED_FIELD = obfrm.NumpyToMEDField_VRij(
                self.input_array[..., 0],
                MED_FIELD,
                BxMesh,
                self.NomCompo,
                temps=[0, 0, 0],
            )
            ml.WriteFieldUsingAlreadyWrittenMesh(fmed_out, MED_FIELD)

        if temps == "instantanne":
            # On va ecrire a chaque pas de temps
            for i, t in enumerate(self.ListeTemps):
                temps = [t, i, 0]
                MED_FIELD = obfrm.NumpyToMEDField_VRij(
                    self.input_array[..., i, 0], MED_FIELD, BxMesh, self.NomCompo, temps
                )
                ml.WriteFieldUsingAlreadyWrittenMesh(fmed_out, MED_FIELD)

    ##### FIN MODIF GAB

    @classmethod
    def initgravity(cls, const, name, other):
        """
        Return a vector Field initialize with three constance

        Args:
            const: a list of three scalar number
            name: a name
            other: a 'model Field' for x axis and labels

        """
        # var = VectorField(name, other._ax, other.labelRef, other.tex)
        c = len(other._ax[0, :])
        var = Field(
            np.zeros((3, c)),
            name=name,
            ax=other._ax,
            label_ref=other.labelRef,
            tex=other.tex,
        )
        var[0, :] = const[0]
        var[1, :] = const[1]
        var[2, :] = const[2]
        return var

    # @classmethod
    # def initsource(cls, const, name, other):
    #     """
    #     Return a scalar field

    #     Args:
    #         const: scalar number
    #         name: a name
    #         other: a 'model Field' for x axis and labels

    #     """
    #     var = VectorField(name, other._ax, other.labelRef, other.tex)
    #     c = len(other._ax[0, :])
    #     var = np.zeros((3, c))
    #     var[0, :] = const[0]
    #     var[1, :] = const[1]
    #     var[2, :] = const[2]
    #     return var

    # TODO: voir si ça sert encore après avoir rajouté l'héritage numpy...
    # def initFromAnalytic(cls, genre, a, b, N, name, label, labelRef, tex):
    #     """
    #     initialize field with analytical expression

    #     Args:
    #         genre: could be S, V, T, T3 or T4 for scalar, vector, Tensor, 3rd order tensor and 4th order tensor
    #         a: =x(0)
    #         b: =x(N)
    #         N: =number of point
    #         name: a name corresponding to a register case
    #         label: a label
    #         labelRef: a label for the x axis
    #         tex: a name for the legend

    #     """
    #     if genre == 'S':
    #         var = ScalarField(name, np.linspace(a, b, N).reshape(1, N), labelRef, tex)
    #         c = len(var._ax[0, :])
    #         if name == 'zero':
    #             var = 0 * np.sin(var._ax)
    #         if name == 'un':
    #             var = np.sin(var._ax) / np.sin(var._ax)
    #         if name == 'sin':
    #             var = np.sin(var._ax)
    #         if name == '2sin':
    #             var = 2 * np.sin(var._ax)
    #         if name == 'sin**2':
    #             var = np.sin(var._ax) * np.sin(var._ax)
    #         if name == 'cos':
    #             var = np.cos(var._ax)
    #         if name == 'tan':
    #             var = np.tan(var._ax)
    #     if genre == 'V':
    #         var = VectorField(name, np.linspace(a, b, N).reshape(1, N), labelRef, tex)
    #         c = len(var._ax[0, :])
    #         var = np.zeros((3, c))
    #         if name == 'poiseuille':
    #             var[0, :] = -0.5 * (2 * var._ax - var._ax * var._ax)
    #         if name == 'poly':
    #             var[0, :] = var._ax[:]
    #             var[1, :] = 0.5 * var._ax * var._ax
    #             var[2, :] = 0.166666 * var._ax * var._ax * var._ax
    #         if name == 'dpoly':
    #             var[0, :] = 1
    #             var[1, :] = var._ax[:]
    #             var[2, :] = 0.5 * var._ax * var._ax
    #         if name == 'titi':
    #             var[0, :] = 1
    #             var[1, :] = np.cos(var._ax)
    #             var[2, :] = np.sqrt(np.abs(var._ax))
    #     if genre == 'T':
    #         var = Tensor2Field(name, np.linspace(a, b, N).reshape(1, N), labelRef, tex)
    #         c = len(var._ax[0, :])
    #         var = np.zeros((3, 3, c))
    #         if name == 'identite':
    #             a = np.zeros((3, c))
    #             a[0, 0, :] = 1.0
    #             a[1, 1, :] = 1.0
    #             a[2, 2, :] = 1.0
    #         if name == 'VV':
    #             a = np.zeros((3, c))
    #             a[0, :] = var._ax[:]
    #             a[1, :] = 0.5 * var._ax * var._ax
    #             a[2, :] = 0.166666 * var._ax * var._ax * var._ax
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, :] = a[i, :] * a[j, :]
    #         if name == 'sin':
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, :] = np.sin(var._ax)
    #         if name == 'const':
    #             var[0, 0, :] = 1
    #             var[1, 0, :] = 2
    #             var[2, 0, :] = 3
    #             var[0, 1, :] = 4
    #             var[1, 1, :] = 5
    #             var[2, 1, :] = 6
    #             var[0, 2, :] = 7
    #             var[1, 2, :] = 8
    #             var[2, 2, :] = 9
    #         if name == 'transpo':
    #             var[0, 0, :] = 1
    #             var[0, 1, :] = 2
    #             var[0, 2, :] = 3
    #             var[1, 0, :] = 4
    #             var[1, 1, :] = 5
    #             var[1, 2, :] = 6
    #             var[2, 0, :] = 7
    #             var[2, 1, :] = 8
    #             var[2, 2, :] = 9
    #         if name == 'dconst':
    #             var[0, 0, :] = 0
    #             var[1, 0, :] = 0
    #             var[2, 0, :] = 0
    #             var[0, 1, :] = 0
    #             var[1, 1, :] = 0
    #             var[2, 1, :] = 0
    #             var[0, 2, :] = 0
    #             var[1, 2, :] = 0
    #             var[2, 2, :] = 0
    #         if name == 'poiseuille':
    #             var[0, 2, :] = -0.5 * (2 - 2 * var._ax)

    #         if name == 'ijkk':
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, :] = (i + 1) + j * 3 + (i + 1) + j * 3 + 9 + 27 + (i + 1) + j * 3 + 18 + 54
    #         if name == 'kijk':
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, :] = 1 + i * 3 + j * 9 + 2 + i * 3 + j * 9 + 27 + 3 + i * 3 + j * 9 + 2 * 27
    #         if name == 'ikjk':
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, :] = 1 + j * 3 + i * 9 + 2 + j * 3 + i * 9 + 27 + 3 + j * 3 + i * 9 + 2 * 27
    #         if name == 'ikjk':
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, :] = (i + 1) + j * 9 + (i + 1) + 3 + j * 9 + 27 + (
    #                             i + 1) + 2 * 3 + j * 9 + 2 * 27

    #     if genre == 'T3':
    #         var = Tensor3Field(name, np.linspace(a, b, N).reshape(1, N), labelRef, tex)
    #         c = len(var._ax[0, :])
    #         var = np.zeros((3, 3, 3, c))
    #         if name == 'gradsin':
    #             for i in dims:
    #                 for j in dims:
    #                     var[i, j, 2, :] = np.cos(var._ax)
    #         if name == 'sin':
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         var[i, j, k, :] = np.sin(var._ax)

    #         if name == 'TV':
    #             a = np.zeros((3, c))
    #             a[0, :] = var._ax[:]
    #             a[1, :] = 0.5 * var._ax * var._ax
    #             a[2, :] = 0.166666 * var._ax * var._ax * var._ax
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         for l in dims:
    #                             var[i, j, k, :] = ((i + 1) + j * 3) * a[k, :]
    #         if name == 'VT':
    #             a = np.zeros((3, c))
    #             a[0, :] = var._ax[:]
    #             a[1, :] = 0.5 * var._ax * var._ax
    #             a[2, :] = 0.166666 * var._ax * var._ax * var._ax
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         for l in dims:
    #                             var[i, j, k, :] = ((j + 1) + k * 3) * a[i, :]

    #         if name == 'file':
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         var[i, j, k, :] = (i + 1) + j * 3 + k * 9

    #     if genre == 'T4':
    #         var = Field(name, np.linspace(a, b, N).reshape(1, N), labelRef, tex)
    #         c = len(var._ax[0, :])
    #         var = np.zeros((3, 3, 3, 3, c))
    #         if name == 'cte':
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         for l in dims:
    #                             var[i, j, k, l, :] = (i + 1) + j * 3 + k * 9 + l * 27

    #         if name == 'gradsin':
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         var[i, j, k, 2, :] = np.cos(var._ax)
    #         if name == 'TT':
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         for l in dims:
    #                             var[i, j, k, l, :] = ((i + 1) + j * 3) * ((k + 1) + l * 3)

    #         if name == 'T3V':
    #             a = np.zeros((3, c))
    #             a[0, :] = var._ax[:]
    #             a[1, :] = 0.5 * var._ax * var._ax
    #             a[2, :] = 0.166666 * var._ax * var._ax * var._ax
    #             for i in dims:
    #                 for j in dims:
    #                     for k in dims:
    #                         for l in dims:
    #                             var[i, j, k, l, :] = np.sin(var._ax) * a[l, :]
    #     return var

    def plot(
        self,
        *args,
        same_ax=None,
        fig_num=None,
        subplots=None,
        show=True,
        index=None,
        delta=1.0,
        add_legend=None,
        style=None,
        legend_type=None,
        **kwargs
    ):
        """
        Plot self in an adaptative way.

        Args:
            same_ax (list of tuples):
            fig_num (int): figure number
            subplots (SubplotAxes): subplots to reuse
            show (bool): either to show immedialy
            index (list of list): list of indexes to use
            delta (float): :math:`\\delta` half width of the channel
            add_legend (str): string to add to every legend
            style (str): style for matplotlib to use
            legend_type (dict): {tex: False, val:False, ...}

        Returns:
            SubplotAxes: the subplot corresponding to the plot

        Examples:

            .. jupyter-execute::

                from commons.Tools import *
                a = Field(np.arange(3*15).reshape((3,15)),
                          ax=np.arange(15).reshape((1,15)), name='Vector')
                %matplotlib inline
                ax = a.plot()

        """
        if legend_type is None:
            legend_type = {"val": False, "name": True, "tex": True}
        if add_legend is None:
            add_legend = ""
        elif isinstance(add_legend, Field):
            add_legend = r", " + add_legend.get_tex(val=True, tex=False)
        if same_ax is None:
            same_ax = []

        # setting index, figsize, subplot_x and subplot_y
        dims = len(self.shape)
        if ((dims == 2) and (self.shape[0] == 1)) or (dims == 1):
            if index is None:
                index = [[0, 0]]
            subplot_x = 1
            subplot_y = 1
            figsize = (9, 4)
        elif dims == 2:
            if index is None:
                index = [[0, i] for i in range(self.shape[0])]
            figsize = (9, 3 * self.shape[0] + 1)
            subplot_x = 1
            subplot_y = self.shape[0]
        elif dims == 3:
            if index is None:
                index = [
                    [i, j] for i in range(self.shape[0]) for j in range(self.shape[1])
                ]
            subplot_x = self.shape[0]
            subplot_y = self.shape[1]
            figsize = (9, 3 * self.shape[1] + 1)
        else:
            raise Exception("Dimension is too high, not implemented yet")

        # printing self.tex
        # self.get_tex(p=True)

        # mapping
        ind = self[..., 0].size
        mapping, n_plot = map_same_ax(same_ax, ind)
        if subplots is None:
            plt.figure(fig_num, figsize=figsize)
            subplots = {}
            for j in range(1, n_plot + 1):
                subplots[j] = plt.subplot(subplot_y, subplot_x, j)

        # plotting
        k = 1
        for i, j in index:
            if dims == 1:
                ind = "$"
                if self.shape[0] == 1:
                    Y = np.array([self._npa] * self._ax.shape[1])
                else:
                    Y = self
            elif (dims == 2) and (self.shape[0] == 1):
                ind = "$"
                Y = self[0]
            elif dims == 2:
                ind = r", ind = %d$" % (j + 1)
                Y = self[j]
            elif dims == 3:
                ind = r", ind = %d, %d$" % (i + 1, j + 1)
                Y = self[i, j]
            else:
                raise Exception("NotImplemented")

            # subplots[mapping[k]].title(r'$%s_{%s}$' % (self.name, ind))
            if len(self.label) != 0:
                label = "$" + self.label[k - 1] + "$"
            else:
                label = self.get_tex(**legend_type)[:-1] + ind
            if style is None:
                subplots[mapping[k]].plot(
                    self._ax[0] / delta, Y, label=label + add_legend, **kwargs
                )
            else:
                subplots[mapping[k]].plot(
                    self._ax[0] / delta, Y, style, label=label + add_legend, **kwargs
                )
            if self.labelRef is not None:
                if isinstance(delta, Field):
                    labelref = (
                        r"$"
                        + self.labelRef
                        + r"/ "
                        + delta.get_tex(tex=False, without_dollar=True)
                        + r"$"
                    )
                elif delta != 1.0:
                    labelref = r"$" + self.labelRef + r"/ \delta$"
                else:
                    labelref = r"$" + self.labelRef + r"/ ? $"
                subplots[mapping[k]].set_xlabel(labelref, size="x-large")
            subplots[mapping[k]].legend()
            subplots[mapping[k]].minorticks_on()
            subplots[mapping[k]].grid(True, which="major")
            subplots[mapping[k]].grid(True, which="minor", alpha=0.2)
            k = k + 1

        plt.tight_layout()

        # showing
        if show:
            plt.show()

        return subplots

    def plot2d(
        self,
        *args,
        same_ax=None,
        fig_num=None,
        subplots=None,
        show=True,
        index=None,
        deltax=1.0,
        deltay=1.0,
        add_legend=None,
        style=None,
        n_levels=10,
        **kwargs
    ):
        """
        Plot self in an adaptative way.

        Args:
            same_ax (list of tuples):
            fig_num (int): figure number
            subplots (SubplotAxes): subplots to reuse
            show (bool): either to show immedialy
            index (list of list): list of indexes to use

        Returns:
            SubplotAxes: the subplot corresponding to the plot

        Examples:

            .. jupyter-execute::

                from commons.Tools import *
                grid_x, grid_y = np.mgrid[0:3, 0:5]
                mesh = np.vstack((grid_x.flatten(), grid_y.flatten()))
                print(mesh)
                a = Field(np.arange(3*15).reshape((3,15)),
                          ax=mesh, name=r'\\textrm{Vector over a 2d mesh}')
                %matplotlib inline
                ax = a.plot2d()

        """
        if add_legend is None:
            add_legend = ""
        elif isinstance(add_legend, Field):
            add_legend = r", " + add_legend.get_tex(val=True, tex=False)
        if same_ax is None:
            same_ax = []

        # setting index, figsize, subplot_x and subplot_y
        dims = len(self.shape)
        if ((dims == 2) and (self.shape[0] == 1)) or (dims == 1):
            if index is None:
                index = [[0, 0]]
            subplot_x = 1
            subplot_y = 1
            figsize = (9, 4)
        elif dims == 2:
            if index is None:
                index = [[0, i] for i in range(self.shape[0])]
            figsize = (9, 3 * self.shape[0] + 1)
            subplot_x = 1
            subplot_y = self.shape[0]
        elif dims == 3:
            if index is None:
                index = [
                    [i, j] for i in range(self.shape[0]) for j in range(self.shape[1])
                ]
            subplot_x = self.shape[0]
            subplot_y = self.shape[1]
            figsize = (9, 3 * self.shape[1] + 1)
        else:
            raise Exception("Dimension is too high, not implemented yet")

        # printing self.tex
        # self.get_tex(p=True)

        # mapping
        ind = self[..., 0].size
        mapping, n_plot = map_same_ax(same_ax, ind)
        if subplots is None:
            plt.figure(fig_num, figsize=figsize)
            subplots = {}
            # TODO: remplacer subplot par subplots, pour utiliser de manière cohérente l'API O-O matplotlib
            for j in range(1, n_plot + 1):
                subplots[j] = plt.subplot(subplot_y, subplot_x, j)

        # plotting
        k = 1
        for i, j in index:
            if dims == 1:
                ind = "$"
                Y = self
            elif (dims == 2) and (self.shape[0] == 1):
                ind = "$"
                Y = self[0]
            elif dims == 2:
                ind = r", ind = %d$" % (j + 1)
                Y = self[j]
            elif dims == 3:
                ind = r", ind = %d, %d$" % (i + 1, j + 1)
                Y = self[i, j]
            else:
                raise Exception("NotImplemented")

            # subplots[mapping[k]].title(r'$%s_{%s}$' % (self.name, ind))
            if len(self.label) != 0:
                label = "$" + self.label[k - 1] + "$"
            else:
                label = self.get_tex()[:-1] + ind
            Y = Y.view(np.ndarray)
            mini = np.min(Y)
            maxi = np.max(Y)
            if style is None:
                clb = subplots[mapping[k]].tricontourf(
                    self._ax[0] / deltax,
                    self._ax[1] / deltay,
                    Y,
                    levels=np.linspace(mini, maxi, n_levels),
                    **kwargs
                )
                plt.colorbar(clb, ax=subplots[mapping[k]])
            else:
                clb = subplots[mapping[k]].tricontourf(
                    self._ax[0] / deltax,
                    self._ax[1] / deltay,
                    Y,
                    style,
                    levels=np.linspace(mini, maxi, n_levels),
                    **kwargs
                )
                plt.colorbar(clb, ax=subplots[mapping[k]])
            subplots[mapping[k]].set_title(label + add_legend)
            if self.labelRef is not None:
                if deltax != 1.0:
                    labelref = r"$" + self.labelRef + r"/ \delta$"
                else:
                    labelref = "$" + self.labelRef + "$"
                subplots[mapping[k]].set_xlabel(labelref, size="x-large")
            k = k + 1

        plt.tight_layout()

        # showing
        if show:
            plt.show()

        return subplots

    def inv(self, tol=0.0):
        """
        Return the inverse of self
        """
        a = len(self._npa[0, :])
        b = len(self._npa[:, 0])
        grad = self.createField(np.zeros((b, a)), name="inv")
        for i in range(a):
            for j in range(b):
                if abs(self._npa[j, i]) < tol:
                    grad._npa[j, i] = 1.0 / tol
                    # raise Exception('Warning inverse of %s below tolerance %g' % (self.name, tol))
                else:
                    grad._npa[j, i] = 1.0 / self._npa[j, i]

        return grad


def map_same_ax(list_of_tuple, ind):
    """
    Map the subplot indices according to ``list_of_tuple``

    Args:
        list_of_tuple (list of tuples): ex: [(1,2),(5,4)]
        ind (int): the number of subplots to map

    Returns:
        dict: the mapping

    Examples:
        >>>map_same_ax([(1,2), (5,4)], 6)
        {1: 1, 2: 1, 3: 2, 5: 3, 4: 3, 6: 4}

    """
    dic = {}
    count = 1
    isintup = [False] * ind
    for k in range(1, ind + 1):
        if not isintup[k - 1]:
            for tup in list_of_tuple:
                if k in tup:
                    for j in tup:
                        dic[j] = count
                        isintup[j - 1] = True
                    count += 1
        if not isintup[k - 1]:
            dic[k] = count
            count += 1
    return dic, count - 1


def fluctuijk(uuu, u, uu):
    """
    return the mean field of the correlation fluctuation of three variables.
    @param mean of AAA
    @param mean of A
    @param mean of AA
    return mean of A'A'A'
    """
    Rijk = uuu
    Rij_uk = u.tensorProduct(uu)
    # print(Rij_uk)
    uijk = u.tensorProduct(u).tensorProduct(u)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # GB : For a true permutation, it should not be Rij_uk._npa[j,i,k,:] but Rij_uk._npa[j,k,i,:]
                # But it doesn't matter because the permutation (i,k) -> (k,i) plays on the Rij part (last 2 indices) that are symmetric.
                Rijk._npa[i, j, k, :] = (
                    Rijk._npa[i, j, k, :]
                    + uijk._npa[i, j, k, :] * 2.0
                    - Rij_uk._npa[j, i, k, :]
                    - Rij_uk._npa[i, j, k, :]
                    - Rij_uk._npa[k, i, j, :]
                )
    return Rijk


def fluctuij(p, du, pdu, I=None):
    """
    return the mean field of the correlation fluctuation of two variables.
    @param mean of A
    @param mean of B
    @param mean of AB
    return mean of A'B'
    """
    if (p.TypeField == "vector") and (du.TypeField == "vector"):
        pdu_bis = p.MatProduct(du)
    elif (p.TypeField == "scalar") and (du.TypeField == "vector"):
        pdu_bis = p * du
    elif (p.TypeField == "scalar") and (du.TypeField == "tensor2"):
        pdu_bis = du * p
    elif (p.TypeField == "tensor2") and (du.TypeField == "tensor2"):
        pdu_bis = p.tensorProduct(du)
        pdu_bis = pdu_bis.contract("ij", "ikjk")
    else:
        raise NotImplementedError
    if I is not None:
        pdu_bis = pdu_bis * (I.inv(0.0001))
    res = pdu - pdu_bis
    return res


def load(fichier):
    """
    load a dictionnary
    @param a file name
    """
    output = open(fichier, "rb")
    Run = pickle.load(output, fix_imports=True)
    output.close()
    clef(Run)
    return Run


def loadSansClef(fichier):
    """
    load a dictionnary without keys
    @param a file name
    """
    output = open(fichier, "rb")
    Run = pickle.load(output)
    output.close()
    return Run


def save(Var, fichier):
    """
    save a dictionnary
    @param a dictionnary
    @param a file name
    """

    output = open(fichier, "wb")
    pickle.dump(Var, output)
    output.close()
    return


def clef(SphericalRun):
    """
    print(all the keys of a dictionnary)
    @param a dictionnary
    """
    print("")
    print("Liste des objets")
    for cle in SphericalRun.keys():
        i = 0
        try:
            for cle2 in SphericalRun[cle].keys():
                i = i + 1
                if i == 1:
                    print(cle + "--->" + cle2)
                else:
                    print("------->" + cle2)
        except:
            print("--->" + cle)

    print("")
    return


def NetContribution(Mk, db=0.3, L=2.0):
    """
    return mean value without taking into account 0.
    @param the Field (vector field)
    """
    dk = Mk._ax[0, 2] - Mk._ax[0, 1]
    MkNet = Mk * 0.0
    for j in range(len(Mk._npa[:, 0])):
        for i in range(len(Mk._npa[j, :])):
            som = 0.0
            somax = 0.0
            #### trouver les bornes
            Ln = len(Mk._npa[j, :])
            dbn = int(Ln * db / L)
            indix = []
            for k in range(int(dbn / 2)):
                if k == 0:
                    indix.append(i - k)
                    continue
                if i - k in range(Ln):
                    indix.append(i - k)
                if i + k in range(Ln):
                    indix.append(i + k)

            indix.sort()
            for k in indix[:-1]:
                som = som + (0.5 * Mk._npa[j, k] + 0.5 * Mk._npa[j, k + 1])
                somax = somax + 1.0

            # print(i, len(indix), len(indix[:-1]), somax, som, Mk._npa[j,i], som+Mk._npa[j,i])
            if somax == 0:
                MkNet._npa[j, i] = 0.0
            else:
                # TODO: is it really somax ?
                MkNet._npa[j, i] = som / somax

            # print(i, somax)
    return MkNet


def absoluteValue(Mk):
    """
    return the absolute value
    @param the Field (vector field)
    """
    MkNet = Mk * 0.0
    for j in range(len(Mk._npa[:, 0])):
        for i in range(len(Mk._npa[j, :])):
            MkNet._npa[j, i] = abs(Mk._npa[j, i])
    return MkNet


def cell_to_cell_gradient_scalar(field, dfield, bc_type=0):
    """
    return the first derivative of field
    @param the field. Component by component
    @param the axis field
    @param bc_type for other boundary condition
    """
    # dfield = np.squeeze(dfield)
    # field = np.squeeze(field)
    gradi = np.zeros(len(field))
    # gradi_t = np.zeros(len(field))
    print("field", field.shape)
    print("dfield", dfield.shape)
    gradi[1:-1] = (field[2:] - field[:-2]) / (dfield[2:] - dfield[:-2])
    # gradi_t[1:-1] = (field[2:] - field[:-2]) / (dfield[2:] - dfield[:-2])

    if bc_type == 1:
        # super bizarre cette condition et cette boucle...
        for i in range(int(len(gradi[1:-1]) / 4)):
            if i == 0:
                gradi[i] = (
                    -0.166666666 * field[i + 6]
                    + 1.2 * field[i + 5]
                    - 3.75 * field[i + 4]
                    + 6.66666666 * field[i + 3]
                    - 7.5 * field[i + 2]
                    + 6 * field[i + 1]
                    - 2.45 * field[i]
                ) / (dfield[i + 1] - dfield[i])
            else:
                gradi[i] = (
                    -0.166666666 * field[i + 6]
                    + 1.2 * field[i + 5]
                    - 3.75 * field[i + 4]
                    + 6.66666666 * field[i + 3]
                    - 7.5 * field[i + 2]
                    + 6 * field[i + 1]
                    - 2.45 * field[i]
                ) / (dfield[i + 1] - dfield[i])
                gradi[-i] = (
                    0.166666666 * field[-i - 6]
                    - 1.2 * field[-i - 5]
                    + 3.75 * field[-i - 4]
                    - 6.66666666 * field[-i - 3]
                    + 7.5 * field[-i - 2]
                    - 6 * field[-i - 1]
                    + 2.45 * field[-i]
                ) / (dfield[-i] - dfield[-i - 1])

        # TODO: comprendre cette CL et faire sauter cette boucle for qui m'a l'air très étrange
        # qu = int(len(gradi[1:-1]) / 4)
        # gradi_t[0] = (-0.166666666 * field[6] + 1.2 * field[5] - 3.75 * field[4] + 6.66666666 * field[3] - 7.5 *
        #             field[2] + 6 * field[1] - 2.45 * field[0]) / (dfield[1] - dfield[0])
        #
        # gradi_t[1:qu] = (-0.166666666 * field[7:qu+6] + 1.2 * field[6:qu+5] - 3.75 * field[5:qu+4] + 6.66666666 * field[
        #         4:qu+3] - 7.5 * field[3:qu+2] + 6 * field[2:qu+1] - 2.45 * field[1:qu]) / (
        #         dfield[2:qu+1] - dfield[1:qu])
        #
        # gradi_t[-qu+1:] = (0.166666666 * field[-i - 6] - 1.2 * field[-i - 5] + 3.75 * field[-i - 4] - 6.66666666 *
        #         field[-i - 3] + 7.5 * field[-i - 2] - 6 * field[-i - 1] + 2.45 * field[-i]) / (
        #         dfield[-i] - dfield[-i - 1])

    elif bc_type == 0:
        gradi[0] = (
            -0.166666666 * field[6]
            + 1.2 * field[5]
            - 3.75 * field[4]
            + 6.66666666 * field[3]
            - 7.5 * field[2]
            + 6 * field[1]
            - 2.45 * field[0]
        ) / (dfield[1] - dfield[0])
        gradi[-1] = (
            0.166666666 * field[-7]
            - 1.2 * field[-6]
            + 3.75 * field[-5]
            - 6.66666666 * field[-4]
            + 7.5 * field[-3]
            - 6 * field[-2]
            + 2.45 * field[-1]
        ) / (dfield[-1] - dfield[-2])

    elif bc_type == 2:  # condition periodique
        gradi[0] = (field[1] - field[-1]) / (dfield[1] - dfield[-1])
        gradi[-1] = (field[0] - field[-2]) / (dfield[0] - dfield[-2])

    return gradi


def addpar(string1, string2, sym):
    if sym == " + ":
        return string1 + sym + addparp(string2)
    elif sym == " - ":
        return string1 + sym + addparf(string2)
    elif sym == " / ":
        return "\\frac{" + string1 + "}{" + string2 + "}"
    elif sym == " * ":
        return addparf(string1) + " " + addparf(string2)
    elif sym == "}^{":
        return "{" + string1 + sym + string2 + "}"
    elif sym == " ? ":
        return "(" + string1 + ")" + sym + "(" + string2 + ")"
    else:
        raise Exception("Problème ce symbole n'est pas reconnu")


def addparf(string):
    """
    Cette fonction ajouter des parentheses s'il y en a besoin pour une opération de type multiplication à la chaine de charactère.

    Args:
        string (str): la chaine a évaluer

    Return:
        str: la chaine modifiée
    """
    add = False
    openpar = 0
    for char in string:
        if ((char == "+") or (char == "-")) and (openpar == 0):
            add = True
        if char == "(" or char == "{":
            openpar += 1
        elif char == ")" or char == "}":
            openpar -= 1
    if add:
        string = "(" + string + ")"
    return string


def addparp(string):
    """
    Cette fonction décide s'il faut ajouter des parenthèses après une opération +
    """
    add = False
    if string[0] == "-":
        add = True
    if add:
        string = "(" + string + ")"
    return string


# je l'ai rajoutée parce que cette méthode est appelée par moentum_liq
def tracer(
    bfield,
    title,
    axe=None,
    textitle="",
    ylim=None,
    xlim=None,
    grid=False,
    legend="outside",
    markerSize=None,
    log=False,
    couleur=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    markevry=[3],
    yscale=1.0,
    xscale=1.0,
    scalepng=1.5,
    applatir=False,
    plotall=True,
    sym=[10],
    zerocentre=False,
    doubleLegend=False,
    casename=None,
):
    """
    Creation of a parametrized plotinitialization of a field object

    Args:
        bfield: the field's list to plot (can be of any sub-field type)
        title: the name of the png file created
        axe: can contain the list of components to be plotted. Default is None for all
        textitle: the title of the plot (optional)
        ylim: limit of y axis. Default is None for all
        xlim: limit of x axis. Default is None for all
        grid: Plot the grid background. Default is True
        legend: ... pas au point
        markerSize: Default is an autoscale from the scale of the png file (scalepng).
        log: Plot with logarithmic scale. Default is False
        couleur: ex couleur=[3,2,3,1] means the first 3 variables are plotted in the same color, then the following two
            with an other etc..
        markevry:
        yscale: dimenssionless tools for y axis
        xscale: dimenssionless tools for x axis
        scalepng: scale of the out.png
        applatir: aspect ratio of the picture = 0.5. Default is False
    """

    if len(markevry) == 1:
        cte = markevry[0]
        tmp = [
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
            cte,
        ]
        markevry = tmp
    print("exporting graph " + title + ".png")
    global compteur
    compteur = compteur + 1
    f = plt.figure(compteur)

    linestyles = []
    markers = [
        "d",
        ".",
        ".",
        ".",
        ".",
        None,
        ".",
        ".",
        ".",
        ",",
        "x",
        ".",
        "1",
        "p",
        "3",
        "2",
        "4",
        "H",
        "v",
        "8",
        "<",
    ]
    markerfacecols = [
        "white",
        "full",
        "full",
        "white",
        "full",
        "full",
        "full",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
        "r",
    ]
    markersizes = [
        2.5,
        4.0,
        2.0,
        4.0,
        2.0,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
        scalepng * 1.5,
    ]

    markerSize = scalepng * 1.5

    if couleur == "standard":
        colors = ["k", "r", "g", "m", "c", "y", "b", "0.8", "0.5"]
    elif couleur == "grey":
        colors = ["1.0", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"]
    elif couleur == "mychoice1":
        colors = ["r", "c", "m", "g", "y", "k", "b"]
    elif couleur == "red":
        colors = ["r", "r", "r", "r", "r", "r", "r"]
    else:
        colors = []
        colori = [
            "k",
            "r",
            "m",
            "g",
            "y",
            "b",
            "0.8",
            "0.5",
        ]  # GB-HERE is default colors
        for k in range(len(couleur)):
            try:
                for i in range(couleur[k]):
                    colors.append(colori[k])
            except:
                colors.append("k")

    colors.append("k")
    styles = markers
    axisNum = 0
    LineStyleNum = 0
    ax = plt.subplot(111)

    if sym == [0]:
        sym_sca = [1]
        sym_vec = [1, 1, 1]
        sym_tens = [1, 1, -1, 1, 1, -1, -1, -1, 1]
    elif sym == [1]:
        sym_sca = [1]
        sym_vec = [1, 1, 1]
        sym_tens = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    elif sym == [-1]:
        sym_sca = [-1]
        sym_vec = [-1, -1, -1]
        sym_tens = [-1, -1, -1, -1, -1, -1, -1, -1, -1]

    #     compt=0
    for i in range(len(bfield)):
        if sym == [10] or sym == 60:
            pass
        else:
            if bfield[i].TypeField == "scalar":
                bfield[i].symetriser(sym_sca, axe)
            elif bfield[i].TypeField == "vector":
                bfield[i].symetriser(sym_vec, axe)
            elif bfield[i].TypeField == "tensor2":
                bfield[i].symetriser(sym_tens, axe)
            else:
                NotImplementedError

        if (bfield[i].TypeField == "tensor2") or (bfield[i].TypeField == "tensor3"):
            bfield[i].preReshape()
        if zerocentre:
            ### mettre zero au centre
            # if isinstance(bfield[i],VectorField):
            #    L=len(bfield[i]._npa[0,:])/2
            #    milieu_x=bfield[i]._npa[0,L]
            #    milieu_y=bfield[i]._npa[1,L]
            #    milieu_z=bfield[i]._npa[2,L]
            #    bfield[i]._npa[0,:]-=milieu_x
            #    bfield[i]._npa[1,:]-=milieu_y
            #    bfield[i]._npa[2,:]-=milieu_z
            # mettre zero en y=0
            if bfield[i].TypeField == "vector":
                milieu_x = bfield[i]._npa[0, 0]
                milieu_y = bfield[i]._npa[1, 0]
                milieu_z = bfield[i]._npa[2, 0]
                bfield[i]._npa[0, :] -= milieu_x
                bfield[i]._npa[1, :] -= milieu_y
                bfield[i]._npa[2, :] -= milieu_z

        # j sont les composantes?
        if axe != None:
            if log == True:
                for j in axe:
                    color = colors[axisNum % len(colors)]
                    if (
                        (colors[axisNum - 1] != colors[axisNum] or axisNum == 0)
                        and doubleLegend
                    ) or doubleLegend == False:
                        labell = bfield[i].tex
                    else:
                        labell = ""
                    if LineStyleNum < len(linestyles):
                        ax.semilogx(
                            bfield[i]._ax[0, :] * xscale,
                            bfield[i]._npa[j, :] * yscale,
                            linestyles[LineStyleNum % len(linestyles)],
                            label=labell,
                            color=color,
                            ms=markerSize,
                            markevery=markevry[i],
                            markerfacecolor=color,
                            markeredgecolor=color,
                        )
                    else:
                        style = styles[(LineStyleNum - len(linestyles)) % len(styles)]
                        markerfacecol = markerfacecols[
                            (LineStyleNum - len(linestyles)) % len(styles)
                        ]
                        if markerfacecol == "full":
                            markerfacecol = color
                        markersizess = markersizes[
                            (LineStyleNum - len(linestyles)) % len(styles)
                        ]
                        ax.semilogx(
                            bfield[i]._ax[0, :] * xscale,
                            bfield[i]._npa[j, :] * yscale,
                            linestyle="-",
                            linewidth=0.4,
                            marker=style,
                            label=labell,
                            color=color,
                            ms=markersizess,
                            markevery=markevry[i],
                            markerfacecolor=markerfacecol,
                            markeredgecolor=color,
                            markeredgewidth=0.2,
                        )

                    if colors[axisNum] != colors[axisNum + 1] and couleur != "standard":
                        axisNum += 1
                        LineStyleNum = 0
                    else:
                        axisNum += 1
                        LineStyleNum += 1
            else:
                for j in axe:
                    color = colors[axisNum % len(colors)]
                    if (
                        (colors[axisNum - 1] != colors[axisNum] or axisNum == 0)
                        and doubleLegend
                    ) or doubleLegend == False:
                        labell = bfield[i].tex
                    else:
                        labell = ""
                    if LineStyleNum < len(linestyles):
                        if doubleLegend == True:
                            ax.plot(
                                bfield[i]._ax[0, :] * xscale,
                                bfield[i]._npa[j, :] * yscale,
                                linestyles[LineStyleNum % len(linestyles)],
                                label="",
                                color=color,
                                ms=markerSize,
                                markevery=markevry[i],
                                markerfacecolor=color,
                                markeredgecolor=color,
                            )
                            plt.plot([], [], "-", color=color, label=labell)

                        else:
                            ax.plot(
                                bfield[i]._ax[0, :] * xscale,
                                bfield[i]._npa[j, :] * yscale,
                                linestyles[LineStyleNum % len(linestyles)],
                                label=labell,
                                color=color,
                                ms=markerSize,
                                markevery=markevry[i],
                                markerfacecolor=color,
                                markeredgecolor=color,
                            )
                    else:
                        style = styles[(LineStyleNum - len(linestyles)) % len(styles)]
                        markerfacecol = markerfacecols[
                            (LineStyleNum - len(linestyles)) % len(styles)
                        ]
                        if markerfacecol == "full":
                            markerfacecol = color
                        markersizess = markersizes[
                            (LineStyleNum - len(linestyles)) % len(styles)
                        ]
                        if doubleLegend == True:
                            ax.plot(
                                bfield[i]._ax[0, :] * xscale,
                                bfield[i]._npa[j, :] * yscale,
                                linestyle="-",
                                linewidth=0.4,
                                marker=style,
                                label="",
                                color=color,
                                ms=markersizess,
                                markevery=markevry[i],
                                markerfacecolor=markerfacecol,
                                markeredgecolor=color,
                                markeredgewidth=0.2,
                            )
                            plt.plot([], [], "-", color=color, label=labell)
                        else:
                            ax.plot(
                                bfield[i]._ax[0, :] * xscale,
                                bfield[i]._npa[j, :] * yscale,
                                linestyle="-",
                                linewidth=0.4,
                                marker=style,
                                label=labell,
                                color=color,
                                ms=markersizess,
                                markevery=markevry[i],
                                markerfacecolor=markerfacecol,
                                markeredgecolor=color,
                                markeredgewidth=0.2,
                            )
                    if colors[axisNum] != colors[axisNum + 1] and couleur != "standard":
                        axisNum += 1
                        LineStyleNum = 0
                    else:
                        axisNum += 1
                        LineStyleNum += 1
                        # ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=bfield[i].tex, color=color)
        else:

            if log == True:
                for j in range(len(bfield[i]._npa[:, 0])):
                    if (
                        (colors[axisNum - 1] != colors[axisNum] or axisNum == 0)
                        and doubleLegend
                    ) or doubleLegend == False:
                        labell = bfield[i].tex
                    else:
                        labell = ""
                    color = colors[axisNum % len(colors)]
                    plt.semilogx(
                        bfield[i]._ax[0, :] * xscale,
                        bfield[i]._npa[j, :] * yscale,
                        linestyles[axisNum % len(linestyles)],
                        label=labell,
                        color=color,
                    )
                    axisNum += 1
            else:
                for j in range(len(bfield[i]._npa[:, 0])):
                    if (
                        (colors[axisNum - 1] != colors[axisNum] or axisNum == 0)
                        and doubleLegend
                    ) or doubleLegend == False:
                        labell = bfield[i].tex
                    else:
                        labell = ""
                    color = colors[axisNum % len(colors)]
                    if LineStyleNum < len(linestyles):
                        plt.plot(
                            bfield[i]._ax[0, :] * xscale,
                            bfield[i]._npa[j, :] * yscale,
                            linestyles[axisNum % len(linestyles)],
                            label=labell,
                            color=color,
                            ms=markerSize,
                            markevery=markevry[i],
                            markerfacecolor="none",
                            markeredgecolor=color,
                        )
                    else:
                        style = styles[(LineStyleNum - len(linestyles)) % len(styles)]
                        markerfacecol = markerfacecols[
                            (LineStyleNum - len(linestyles)) % len(styles)
                        ]
                        if markerfacecol == "full":
                            markerfacecol = color
                        markersizess = markersizes[
                            (LineStyleNum - len(linestyles)) % len(styles)
                        ]
                        plt.plot(
                            bfield[i]._ax[0, :] * xscale,
                            bfield[i]._npa[j, :] * yscale,
                            linestyle="-",
                            linewidth=0.4,
                            marker=style,
                            label=labell,
                            color=color,
                            ms=markersizess,
                            markevery=markevry[i],
                            markerfacecolor=markerfacecol,
                            markeredgecolor=color,
                            markeredgewidth=0.2,
                        )

                    if colors[axisNum] != colors[axisNum + 1] and couleur != "standard":
                        axisNum += 1
                        LineStyleNum = 0
                    else:
                        axisNum += 1
                        LineStyleNum += 1

    if doubleLegend:
        if casename == None:
            print("il faut remplir casename si option doublelegend est activee")
        linestyles = ["-"]
        for i in range(len(casename)):
            color = colors[i]
            style = styles[i]
            markerfacecol = markerfacecols[i]
            if markerfacecol == "full":
                markerfacecol = color
            markersizess = markersizes[i]
            if i == 0:
                plt.plot(
                    [],
                    [],
                    linestyles[i],
                    label=casename[i],
                    color="k",
                    linewidth=0.4,
                    marker=style,
                    ms=markersizess,
                    markerfacecolor=markerfacecol,
                    markeredgecolor="k",
                    markeredgewidth=1.0,
                )
            else:
                markerlocal = markers[i - 1]
                if markers[i - 1] == None:
                    markerlocal = "."
                plt.plot(
                    [],
                    [],
                    markerlocal,
                    label=casename[i],
                    color="k",
                    marker=style,
                    ms=markersizess,
                    markerfacecolor=markerfacecol,
                    markeredgecolor="k",
                    markeredgewidth=1.0,
                )

    if grid:
        plt.grid("on")
    if xlim != None:
        ax.set_xlim(xlim, emit=False)
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    else:
        try:
            autoscale_y(ax)
        except:
            print("auto scale fail")
    # plt.xlabel(r'$y/h$')
    ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    plt.ticklabel_format(style="plain", axis="y", scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    dpii = 400
    # plt.title(textitle)
    if applatir:
        f.set_size_inches(scalepng, scalepng / 2.0)
    else:
        f.set_size_inches(scalepng, scalepng)
    inSizeLegend = int(scalepng * 3.5)
    sauv = title + "_off.png"
    if plotall:
        plt.savefig(sauv, bbox_inches="tight", dpi=dpii)
        # plt.legend(loc=(0.2,0.7),prop={'size':inSizeLegend}, frameon=True, ncol=2)
        plt.legend(loc=0, prop={"size": inSizeLegend}, frameon=True, ncol=1)
        sauv = title + "_in.png"
        plt.savefig(sauv, bbox_inches="tight", dpi=dpii)

    plt.legend(loc=0)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc="center left", bbox_to_anchor=(1.1, 0.5), frameon=True)
    sauv = title + "_out.png"
    plt.savefig(sauv, bbox_inches="tight", dpi=dpii)

    for i in range(len(bfield)):
        if (bfield[i].TypeField == "tensor2") or (bfield[i].TypeField == "tensor3"):
            bfield[i].postReshape()
    # plt.close(f)


def tracerJFM(field, name="blabla"):
    print("tracer JFM pictures %s" % name)
    color2 = "#8dd3c7"
    color4 = "#ffffb3"
    color1 = "#bebada"
    color3 = "#fb8072"
    ###### figure 1
    scalepng = 1.5
    f, (ax) = plt.subplots()
    markeredgewidths = [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    ]
    ax.set_xlim([0.0, 1.0], emit=False)

    if "JFMequilibreDNS" in name:
        style = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1000.0),
            (1.0, 1.0),
            (1000.0, 1.0),
            (1.0, 1.0),
            (2.0, 2.0),
        ]
        color = ["k", "g", "r", "b", "m", "c", "k", "g", "r", "b"]
        marker = [None, "d", ".", "x", "d", "+", None, "_", None, None]
        label = [
            r"$\mathbf{\dot M_l=\Sigma M_l}$",
            r"$\mathbf{\dot M_l}$ DNS",
            r"$\mathbf{M_l^\sigma}+\mathbf{M_l^P}$",
            r"$\mathbf{M_l^\Pi}$",
            r"$\mathbf{M_l^\Re}$",
            r"$\mathbf{M_l^\tau}$",
            r"$\alpha_vu_r\wedge(\nabla\wedge u_l)$",
            r"$-\alpha_vu_r\wedge(\nabla\wedge u_l)$",
        ]
        if "JFMequilibreDNSx" in name:
            label = [
                r"$\mathbf{\dot M_l}$",
                r"$\mathbf{M_l^{LD}}$",
                r"$\mathbf{M_l^D}$",
                r"$\mathbf{M_l^{TD}}$",
                r"$\mathbf{M_l^\tau}$",
                r"$\mathbf{M_l^{extra}}$",
                r"$\alpha_v|u_r|u_r$",
                r"$\alpha_v|u_r|u_r$",
            ]

        if "JFMequilibreDNSz" in name:
            label = [
                r"$\mathbf{M_l^{RANS}-\frac{\alpha_l}{\alpha_v}M^L_{\nabla P}}$",
                r"$\frac{\alpha_l}{\alpha_v}\mathbf{M^{LD}}$",
                r"$-\mathbf{M_\Re^L}$",
                r"$\mathbf{-M^{TD}}$",
                r"$\mathbf{-M^\tau}$",
                r"$\mathbf{-M^{extra}}$",
                r"$-\alpha_l\nabla\overline{p_l^{SP}}-\frac{\alpha_l}{\alpha_v}\mathbf{M^L_{\nabla P}}$",
                r"$\alpha_v|u_r|u_r$",
            ]

        if "JFMequilibreDNSx_vapeur" in name:
            label = [
                r"$\mathbf{M_v^{RANS}}-\alpha_v\nabla\left(\overline{p_l^b}-\widetilde{p_l}\right)$",
                r"$\mathbf{M_v^{LD}}$",
                r"$\mathbf{M_v^D}$",
                r"$\mathbf{M_v^{TD}}$",
                r"$\mathbf{M_v^\tau}$",
                r"$\mathbf{M_v^{extra}}$",
                r"$\alpha_v|u_r|u_r$",
                r"$\alpha_v|u_r|u_r$",
            ]

        if "JFMequilibreDNSz_vapeur" in name:
            label = [
                r"$\mathbf{M_v^{RANS}-M^L_{\nabla P}}$",
                r"$\mathbf{M^{LD}}$",
                r"$\mathbf{M^L_\Re}$",
                r"$\mathbf{M^{TD}}$",
                r"$\mathbf{M^\tau}$",
                r"$\mathbf{M^{extra}}$",
                r"$\alpha_v|u_r|u_r$",
                r"$\alpha_v|u_r|u_r$",
            ]

        markerfacecols = [
            "None",
            "none",
            "None",
            "None",
            "None",
            "None",
            "None",
            "None",
            "None",
            "None",
        ]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 3.5,
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
            scalepng * 3,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]

        if "JFMequilibreDNSx" in name:
            ax.set_ylim([-0.015, 0.04], emit=False)  # S180
            # ax.set_ylim([-0.01,0.03], emit=False) #S180g8
            # ax.set_ylim([-0.002,0.01], emit=False) #D180
            # ax.set_ylim([-0.004,0.03], emit=False) #D180g8
            # ax.set_ylim([-0.001,0.004], emit=False) #D127

        if "JFMequilibreDNSz" in name:
            # ax.set_ylim([-0.06,0.06], emit=False) #S180
            # ax.set_ylim([-0.02,0.1], emit=False) #S180g8
            ax.set_ylim([-0.0015, 0.008], emit=False)  # D127
            # ax.set_ylim([-0.002,0.004], emit=False) #D180
            # ax.set_ylim([-0.005,0.02], emit=False) #D180g8

        if "JFMequilibreDNSx_vapeur" in name:
            # ax.set_ylim([-0.01,0.01], emit=False)
            ax.set_ylim([-0.02, 0.015], emit=False)  # S180
            # ax.set_ylim([-0.03,0.01], emit=False) #S180g8
            # ax.set_ylim([-0.01,0.002], emit=False) #D180
            # ax.set_ylim([-0.03,0.003], emit=False) #D180g8
            # ax.set_ylim([-0.004,0.001], emit=False) #D127

        if "JFMequilibreDNSz_vapeur" in name:
            # ax.set_ylim([-0.0003,0.0005], emit=False) #D180
            # ax.set_ylim([-0.0006,0.001], emit=False) #D180g8
            # ax.set_ylim([-0.003,0.003], emit=False) #S180
            # ax.set_ylim([-0.002,0.003], emit=False) #S180g8
            ax.set_ylim([-0.0001, 0.00006], emit=False)  # D127

        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

    if "JFMequilibre_force" in name:
        style = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1.0, 1.0),
            (2.0, 2.0),
        ]
        if "Bois" in name:
            dash = [
                (1000.0, 1.0),
                (1000.0, 1.0),
                (1000.0, 1.0),
                (1000.0, 1.0),
                (1000.0, 1.0),
                (1.0, 1.0),
                (2.0, 2.0),
            ]
        color = ["k", "g", "r", "b", "k", "g", "r", "b", "k", "g"]
        if "Bois" in name:
            color = ["k", "g", "r", "b", "k", "b", "r"]
        marker = ["d", ".", "x", "*", "d", ".", "*", None, None, None]
        if "Bois" in name:
            marker = ["d", ".", "x", "*", "d", None, None, None]

        if "z" in name:
            label = [
                r"$\mathbf{M^{LD}}$",
                r"$\mathbf{M^{TD}}$",
                r"$\mathbf{M^\tau}$",
                r"$\mathbf{M^\Pi}$",
                r"residue",
                r"$-\alpha_v\nabla{\overline{p_l^{SP}}}+\mathbf{M^L_{\nabla P}}$",
                r"$\mathbf{M_\Re^{L}|_{lam}}+\mathbf{M_\Re^{L}|_{turb}}$",
                r"$C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$",
                r"$-C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$",
            ]
        if "Bois" in name:
            label = [
                r"$\mathbf{M^{LD}}$",
                r"$\mathbf{M^{TD}}$",
                r"$-\alpha_v\nabla{\overline{p_l^{SP}}}+\mathbf{M^L_{\nabla P}}$",
                r"$\mathbf{M_\Re^{L}|_{lam}}+\mathbf{M_\Re^{L}|_{turb}}$",
                r"residue",
                r"$C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$",
                r"$-C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$",
            ]

        if "x" in name:
            label = [
                r"$\mathbf{M^{LD}}$",
                r"$\mathbf{M^{TD}}$",
                r"$\mathbf{M^\tau}$",
                r"$\mathbf{M^\Pi}$",
                r"residue",
                r"$-\alpha_v\nabla{\overline{p_l^{,}}}^l$",
                r"$\mathbf{M^{D}}$",
                r"$C_D^*\frac{3}{4d_b}\rho_l\alpha_v|u_r|u_r$",
            ]
        markerfacecols = ["None", "None", "None", "None", "k", "g", "r", "b", "k", "g"]
        markersizes = [
            scalepng * 2,
            scalepng * 3.5,
            scalepng * 2.5,
            scalepng * 3.5,
            scalepng * 2,
            scalepng * 3.5,
            scalepng * 3.0,
            scalepng * 3,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]

        if "x" in name:
            # ax.set_ylim([-0.01,0.01], emit=False)
            # ax.set_ylim([-0.015,0.015], emit=False) #S180
            ax.set_ylim([-0.04, 0.04], emit=False)  # S180g8
            # ax.set_ylim([-0.006,0.006], emit=False) #D180
            # ax.set_ylim([-0.035,0.035], emit=False) #D180g8
            # ax.set_ylim([-0.005,0.005], emit=False) #D127

        if "z" in name:
            # ax.set_ylim([-0.06,0.06], emit=False) #S180
            # ax.set_ylim([-0.02,0.1], emit=False) #S180g8
            # ax.set_ylim([-0.0005,0.003], emit=False) #D127
            ax.set_ylim([-0.002, 0.004], emit=False)  # D180
            # ax.set_ylim([-0.005,0.02], emit=False) #D180g8
            if "vapeur" in name:
                # ax.set_ylim([-0.0003,0.0003], emit=False) #D180
                # ax.set_ylim([-0.0006,0.0006], emit=False) #D180g8
                # ax.set_ylim([-0.003,0.003], emit=False) #S180
                ax.set_ylim([-0.002, 0.002], emit=False)  # S180g8
                # ax.set_ylim([-0.0001,0.0001], emit=False) #D127

        if "lift" in name:
            ax.set_ylim([-0.1, 0.1], emit=False)  # D127

        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

    if "JFMequilibre_force_z_vapeur_2" in name:
        style = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1.0, 1.0),
            (2.0, 2.0),
            (3.0, 1.0),
        ]
        color = ["k", "g", "r", "b", "k", "g", "r", "b", "k", "g"]
        marker = ["d", ".", "x", "*", "d", ".", "*", None, None, None]

        label = [
            r"interface+inter pressure",
            r"reste pressure",
            r"turb disp",
            r"gradP",
            r"lift",
            r"residue",
            r"maintien",
        ]
        markerfacecols = ["None", "None", "None", "None", "k", "g", "r", "b", "k", "g"]
        markersizes = [
            scalepng * 2,
            scalepng * 3.5,
            scalepng * 2.5,
            scalepng * 3.5,
            scalepng * 2,
            scalepng * 3.5,
            scalepng * 3.5,
            scalepng * 3,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_ylim([-70.0, 50.0], emit=False)  # D127
        if "lift" in name:
            ax.set_ylim([-0.1, 0.1], emit=False)  # D127

        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

    if "C_p" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
        ]
        color = ["k", "r", "b", "g", "m", "c"]
        marker = [None, ".", "x", "d", "*", "+"]
        label = [r"DNS", r"$C=1$", r"$C=0.75$", r"$C=0.5$", r"$C=0.25$", r"$C=0$"]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 3,
            scalepng * 2,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]

        if "l" in name:
            # ax.set_ylim([-0.03,0.01], emit=False) #S180g8
            # ax.set_ylim([-0.03,0.01], emit=False) #S180
            # ax.set_ylim([-0.035,0.01], emit=False) #D80g8
            ax.set_ylim([-0.003, 0.001], emit=False)  # D127
            # ax.set_ylim([-0.005,0.0015], emit=False) #D180
        if "v" in name:
            # ax.set_ylim([0.0,3.], emit=False) #S80g8
            # ax.set_ylim([0.0,0.5], emit=False) #S80
            # ax.set_ylim([0.0,0.5], emit=False) #D80g8
            ax.set_ylim([0.0, 0.07], emit=False)  # D127
            # ax.set_ylim([0.0,0.07], emit=False) #D180
        # ax.set_ylim([-0.001,0.001], emit=False) #S180
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

    if "JFMequilibreDNSzapriori" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
        ]
        color = ["k", "b", "r", "k", "b", "r"]
        marker = [None, None, None, "d", "d", "d"]
        label = [
            r"$\mathbf{M^P}$ DNS",
            r"$\mathbf{M^\sigma}$ DNS",
            r"$\mathbf{M^{LD}}$ DNS",
            r"$\mathbf{M^P}$ model",
            r"$\mathbf{M^\sigma}$ model",
            r"$\mathbf{M^{LD}}$ model",
        ]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 2.5,
            scalepng * 2.5,
            scalepng * 2.5,
            scalepng * 2.5,
            scalepng * 2.5,
            scalepng * 2.5,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        # ax.set_ylim([-0.005,0.006], emit=False) #D127
        # ax.set_ylim([-0.015,0.015], emit=False) #D180
        # ax.set_ylim([-0.04,0.05], emit=False) #D180g8
        ax.set_ylim([-0.6, 0.8], emit=False)  # S180g8
        # ax.set_ylim([-0.6,0.8], emit=False) #S180
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

    if "JFMI" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (3.0, 3.0),
        ]
        color = ["k", "b", "r", "g", "m", "k"]
        marker = ["None", "d", ".", "x", "*", None]
        label = [r"D127", r"D180", r"D180g8", r"S180", r"S180g8", ""]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_ylim([0.0, 0.15], emit=False)
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
        ax.set_ylabel(r"$\alpha_v$", fontsize=8)

    if "JFM_I_nept" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [(1000.0, 1.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (1.0, 3.0)]
        color = ["k", "b", "r", "g", "m", "k"]
        marker = ["None", "d", ".", "x", "*", None]
        label = [
            r"DNS",
            r"$Eo_c=10$ standard",
            r"$Eo_c=2.5$ standard",
            r"$Eo_c=10$ model",
            r"$Eo_c=2.5$ model",
        ]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_ylim([0.0, 0.13], emit=False)
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
        ax.set_ylabel(r"$\alpha_v$", fontsize=8)

    if "JFM_I_nept_sphe" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [(1000.0, 1.0), (1.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
        color = ["k", "b", "g", "m", "k"]
        marker = ["None", "d", "x", "*", None]
        label = [r"DNS", r"$Eo_c=10$ standard", r"$Eo_c=10$ model"]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_ylim([0.0, 0.08], emit=False)
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
        ax.set_ylabel(r"$\alpha_v$", fontsize=8)

    if "JFMWe" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (3.0, 3.0),
        ]
        color = ["k", "b", "r", "g", "m", "k"]
        marker = ["None", "d", ".", "x", "*", None]
        label = [r"D127", r"D180", r"D180g8", r"S180", r"S180g8", r"$\mathbf{We_c=3}$"]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_ylim([0.0, 13.0], emit=False)
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
        ax.set_ylabel(r"$\mathbf{We}$", fontsize=8)

    if "JFMvitesse" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (3.0, 3.0),
        ]
        color = ["k", "b", "r", "g", "m", "k"]
        marker = ["None", "d", ".", "x", "*", None]
        label = [r"D127", r"D180", r"D180g8", r"S180", r"S180g8"]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_ylim([-0.1, 1.2], emit=False)
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
        if "_l" in name:
            ax.set_ylabel(r"$\mathbf{\overline{U_l}^l}$", fontsize=8)
        if "_r" in name:
            ax.set_ylabel(r"$\mathbf{\overline{U_v}^v-\overline{U_l}^l}$", fontsize=8)

    if "Rij" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (1000.0, 1.0),
            (3.0, 3.0),
        ]
        color = ["k", "b", "r", "g", "m", "k"]
        marker = ["None", "d", ".", "x", "*", None]
        label = [r"D127", r"D180", r"D180g8", r"S180", r"S180g8"]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
            scalepng * 3,
            scalepng * 3,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
        if "11" in name:
            ax.set_ylabel(r"$\mathbf{R_{11}}$", fontsize=8)
            ax.set_ylim([-0.01, 0.1], emit=False)
            ax.set_xlim([0.0, 2.0], emit=False)
        elif "22" in name:
            ax.set_ylabel(r"$\mathbf{R_{22}}$", fontsize=8)
            ax.set_ylim([0.0, 0.015], emit=False)
        elif "33" in name:
            ax.set_ylabel(r"$\mathbf{R_{22}}$", fontsize=8)
            ax.set_ylim([0.0, 0.015], emit=False)
        elif "13" in name:
            ax.set_ylabel(r"$\mathbf{-R_{12}}$", fontsize=8)
            ax.set_ylim([0.0, 0.007], emit=False)
        else:
            NotImplementedError

    if "Cp" in name:
        style = ["-", "-", "-", "-", "-", "-", "-"]
        dash = [
            (1000.0, 1.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
        ]
        color = ["k", "r", "b", "g", "m", "c"]
        marker = [None, ".", "x", "d", "*", "+"]
        label = [r"DNS", r"$C=1$", r"$C=0.75$", r"$C=0.5$", r"$C=0.25$", r"$C=0$"]
        markerfacecols = ["None", "None", "None", "None", "None", "None"]
        markersizes = [
            scalepng * 3,
            scalepng * 3,
            scalepng * 2,
            scalepng * 2,
            scalepng * 3,
            scalepng * 2,
        ]
        markeredgewidths = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        # ax.set_ylim([-0.005,0.01], emit=False) #D127
        # ax.set_ylim([-0.015,0.03], emit=False) #D180
        # ax.set_ylim([-0.04,0.085], emit=False) #D180g8
        # ax.set_ylim([-0.6,0.9], emit=False) #S180g8
        ax.set_ylim([-0.6, 0.9], emit=False)  # S180
        ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

        if "l" in name:
            # ax.set_ylim([-0.04,0.0], emit=False) #S180g8
            ax.set_ylim([-0.015, 0.01], emit=False)  # S180
            # ax.set_ylim([-0.035,0.01], emit=False) #D80g8
            # ax.set_ylim([-0.004,0.0005], emit=False) #D127
            # ax.set_ylim([-0.005,0.0015], emit=False) #D180
        if "v" in name:
            # ax.set_ylim([0.0,3.], emit=False) #S80g8
            # ax.set_ylim([0.0,0.5], emit=False) #S80
            # ax.set_ylim([0.0,0.4], emit=False) #D80g8
            ax.set_ylim([0.0, 0.04], emit=False)  # D127
            # ax.set_ylim([0.0,0.04], emit=False) #D180

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    # plt.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    plt.tight_layout()
    dpii = 400

    for i in range(len(field)):
        ### avec legende
        # print(np.shape(field[i]._ax[0,:]))
        ax.plot(
            field[i]._ax[0, :],
            field[i]._npa[0, :],
            linestyle=style[i],
            dashes=dash[i],
            linewidth=0.4,
            marker=marker[i],
            label=label[i],
            color=color[i],
            ms=markersizes[i],
            markevery=0.06,
            markerfacecolor=markerfacecols[i],
            markeredgecolor=color[i],
            markeredgewidth=markeredgewidths[i],
        )
        ### sans legende
        # ax.plot(field[i]._ax[0,:], field[i]._npa[0,:], linestyle=style[i], dashes=dash[i], linewidth=0.4, marker=marker[i], color=color[i], ms=markersizes[i], markevery=0.06, markerfacecolor=markerfacecols[i], markeredgecolor = color[i], markeredgewidth=1.0)

    f.set_size_inches(scalepng, scalepng)
    # f.set_size_inches(scalepng*1.5,scalepng) # faire des images rectangulaires
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    ax.legend(loc="center left", bbox_to_anchor=(1.1, 0.5), frameon=True)
    inSizeLegend = int(scalepng * 3.5)
    # plt.legend(loc='center right',prop={'size':inSizeLegend}, frameon=False)
    # plt.legend(loc=(0.5,0.15),prop={'size':inSizeLegend}, frameon=True)
    plt.legend(loc=0, prop={"size": inSizeLegend}, frameon=True, ncol=1)
    # plt.legend().set_visible(False) ### deactive la legende
    # plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True, ncol=2)
    # plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True)
    plt.savefig(name, bbox_inches="tight", dpi=dpii)
    return


def autoscale_y(ax, margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo, hi = ax.get_xlim()
        y_displayed = yd[((xd > lo) & (xd < hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed) - margin * h
        top = np.max(y_displayed) + margin * h
        return bot, top

    lines = ax.get_lines()
    bot, top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot:
            bot = new_bot
        if new_top > top:
            top = new_top

    ax.set_ylim(bot, top)
    return ax


class FictiveBubble(object):
    def __init__(
        self,
        amplitude=0.0,
        istart=0,
        Ist=0.0,
        dIst=0.0,
        ifin=0,
        If=0.0,
        dIf=0.0,
        ipic=0,
        Ip=0.0,
        dIp=0.0,
        alpha=0.0,
        fax=0.0,
        iax=0.0,
        pax=0.0,
    ):
        self.Start = [istart, Ist, dIst, iax]
        self.flatStart = [istart, Ist, dIst, iax]
        self.flatEnd = [istart, Ist, dIst, iax]
        self.fin = [ifin, If, dIf, fax]
        self.pic = [ipic, Ip, dIp, pax]
        self.alpha = alpha
        self._npa = np.array([])
        self.label = np.array([])
        self.amplitude = amplitude

    def findPic(self, I):
        maxi = max(abs(self._npa[:] - 1))
        for i in range(len(self._npa[:])):
            if abs(self._npa[i] - 1) == maxi:
                self.pic[0] = i + self.Start[0]
                self.pic[1] = self._npa[i]
                self.pic[2] = I.grad()._npa[0, i + self.Start[0]]
                self.pic[3] = I._ax[0, i + self.Start[0]]

        print(
            "pic de taux de vide =",
            1 - self.pic[1],
            "soit x=",
            self.pic[3],
            "et dI/dz=",
            self.pic[2],
        )


if __name__ == "__main__":
    c = Field(np.arange(45).reshape((3, 3, 5)), ax=np.arange(5).reshape((1, 5)))
    grad = np.gradient(c, prints=True)
    c.print(field=True)
