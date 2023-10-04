The files *spdata.dat*, *pddata.dat*, and *dfdata.dat* contain the
information in the tables by Antee & O'Mara 1995
(`1995MNRAS.276..859A <https://academic.oup.com/mnras/article/276/3/859/1034639>`_)
, Barklem & O'Mara 1997
(`1997MNRAS.290..102B <https://academic.oup.com/mnras/article/290/1/102/1107836>`_)
, and Barklem et al. 1998
(`1998MNRAS.296.1057B <https://academic.oup.com/mnras/article/296/4/1057/1064804>`_).

Each of these files contain two tables, one for the width cross section
and other for the velocity parameter, one after the other, with the
principal quantum number for the first level (*s*, *p*, and *d* for the
three files, respectively) being the slow rolling index and the
principal quantum number for the second level (*p*, *d*, and *f* for
the three files, respectively) the fast rolling index. The dimenions
of the tables are :math:`21\times18`, :math:`18\times18`, and
:math:`18\times18`,
respectively. The principal quantum numbers and the dimensions are
hardcoded in the *rBarklem* subroutine in the *rbarklem_mod.f90* module
file in the :ref:`Source folder <folder_tree>`. Therefore, in order to
change other than the specific values in the tables, a binary file must
be specified in the :ref:`control file <inputs_control>`.
