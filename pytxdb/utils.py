from sqlalchemy import Table, Column, Float, ForeignKey, Boolean, Date, \
    Integer, String, JSON, create_engine

import pandas as pd
from functools import reduce

from pyranges.pyranges import PyRanges, fill_kwargs, pyrange_apply_single



def check_results(res):
    if len(res) == 0:
        raise ValueError("Did not find any features with those names, "
                         "please check the names you provided")
    else:
        return res



def str_to_type(st, *args):
    """
    return sqlalchemy types based on a string, this is used to create tables dynamically
    :param st: a string
    :param args tablename and column name if you are creating a foreing key
    :return: a sqlachemy type supported types are str, int, float, date, bool and json there is
    limited support for json if not in the list return an error. you can create a foreign key if you provide the
    tablename and column name in *args in that order
    """
    if st == "int":
        coltype = Integer
    elif st == "str":
        coltype = String
    elif st == "float":
        coltype = Float
    elif st == "date":
        coltype = Date
    elif st == "bool":
        coltype = Boolean
    elif st == "json":
        coltype = JSON
    elif st == "fk":
        coltype = ForeignKey("{}.{}".format(args[0], args[1]))
    else:
        raise NotImplementedError("{} is not an implemented column type".format(st))
    return coltype


def dict_to_table(dict, tablename, meta):
    """
    create a sqlalchemy table from a dictionary, currently only a handful of things are supported
    The dictionanry structure is a as follows:
    {colname:
        type: str, int, bool, json, date (only one)
        pk: true
        fk:
            table:table for the foreign key
            column: column name for the foreignkey
    }
    all error checks (like if the table exists etc) are done by sqlalchemy
    :param dict:  a dictionary (see above)
    :param tablename: name of the table
    :param meta: table metadata so they can be added but will not be created until called table.create()
    :return: a sqlalchemy table with metadata (engine data) associated with it
    """
    colnames = list(dict.keys())
    coltypes = []
    for col in colnames:
        if dict[col]["type"] != "fk":
            coltypes.append(str_to_type(dict[col]["type"]))
        else:
            coltypes.append(str_to_type(dict[col]["type"], dict[col]["fk"]["table"], dict[col]["fk"]["column"]))
    idxs = [dict[colname]["index"] for colname in colnames]
    pks = [True if "pk" in dict[colname].keys() else False for colname in colnames]
    table = Table(tablename, meta,
                  *(Column(colname, coltype, primary_key=pk, index=idx)
                    for colname, coltype, pk, idx in zip(colnames, coltypes, pks, idxs))
                  , extend_existing=True)

    return table


def dict_to_engine(params, username=None, pwd=None, host=None, port=None):
    """
    create a database connection based on the yaml description
    :param dict: the output section of the config yaml
    :param kwargs: if type is not sqlite in this order username, password, host, port
    :return: a sqlalchemy engine
    """
    dbtype=params["type"]
    if dbtype == "sqlite":
        dbstring = "sqlite:///{}".format(params["name"])
    elif dbtype == "postgresql":
        dbstring = "postgresql://{}:{}@{}:{}/{}".format(username, pwd, host, port, params["name"])
    elif dbtype == "mariadb":
        dbstring = "mariadb+mariadbconnector:://{}:{}@{}:{}/{}".format(username, pwd, host, port, params["name"])
    elif dbtype == "mysql":
        dbstring = "mysql+mysqlconnector://{}:{}@{}:{}/{}".format(username, pwd, host, port, params["name"])
    else:
        raise NotImplementedError("There is no current support for {}".format(dbtype))
    engine = create_engine(dbstring)
    return engine


def mart_download(mart, fields, common, filters, error="ignore"):
    """
    given a list of annotation types download them from ensembl one by one
    :param mart: biomart dataset connection
    :param fields: fields to download
    :param common: common field like ensembl_gene_id to download with every request and merge
    :param error: what to do if there is an error, "ignore" and continue with the next or
    throw an "error"
    :param: filters: a dict showing which filters to apply this is to prevent new annotations that
    are not in the gtf file being downloaded as well.
    :return: a dataframe of annotations
    """
    annots = []
    for annot in fields:
        if annot == common:
            continue
        else:
            cols = [common, annot]
            try:
                res = mart.query(attributes=[common, annot],
                                 filters=filters)
                res.columns = cols
                annots.append(res)
            except:
                if error=="ignore":
                    print("Could not download {} but will continue with the rest".format(annot))
                else:
                    raise ConnectionError("Could not download {} quitting".format(annot))

    if len(annots)>1:
        annots_merged=reduce(lambda left, right:  # Merge DataFrames in list
               pd.merge(left, right,
                        on=[common],
                        how="left"),
               annots)
    elif len(annots)==1:
        annots_merged=annots.copy()
    else:
        raise ValueError("Could not download any annotations please check your connection")

    return annots_merged

def get_methods(class_instance):
    method_list = [attribute for attribute in dir(class_instance) if callable(getattr(class_instance, attribute))
                   and attribute.startswith('__') is False]
    return method_list

#TODO test some of the used methods in the class methods, expecially the print methods since all the
# genome functions will return a gr object
class GRanges(PyRanges):
    def __init__(self, df=None, chromosomes=None,
                 starts=None, ends=None, strands=None):
        super().__init__(df, chromosomes=chromosomes,
                         starts=starts, ends=ends,
                         strands=strands)

        methods=get_methods(super())
        for method_name in methods:
            method = getattr(super(), method_name)
            self.__setattr__(method_name, classmethod(method))

    def slice(self, start, end):
        """
        row slice ala dataframe
        :param start: start location (not in genomic location but in row location)
        :type start: int
        :param end: end location
        :type end: int
        :param in_place: whether to modify in place
        :type in_place: bool
        :return: sliced version of the GRanges instance
        :rtype: GRanges
        """
        df=self.as_df()
        df=df.iloc[start:end, :]
        new_gr=GRanges(df=df)
        return new_gr








