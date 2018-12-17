[1mdiff --git a/eqfm.py b/eqfm.py[m
[1mindex 87ef08f..bd3f720 100644[m
[1m--- a/eqfm.py[m
[1m+++ b/eqfm.py[m
[36m@@ -53,22 +53,10 @@[m [mclass stateObject(object):[m
         elif type(data) is NoneType:[m
             self.__data = {}[m
     [m
[31m-    def __setitem__(self,key,value):[m
[31m-        check_keys(key,self.__keys)[m
[31m-        self.__data[key] = value[m
[31m-    [m
[31m-    def __getitem__(self,key):[m
[31m-        return self.__data[key][m
[31m-    [m
     @property[m
     def variableDict(self):[m
         return self.__data[m
     [m
[31m-    @variableDict.setter[m
[31m-    def variableDict(self,dataDict):[m
[31m-        assert type(dataDict) is dict[m
[31m-        self.__data = dataDict[m
[31m-    [m
     @property[m
     def category(self):[m
         return self.__category[m
[36m@@ -77,10 +65,6 @@[m [mclass stateObject(object):[m
     def purpose(self):[m
         return self.__purpose[m
     [m
[31m-    def keys(self):[m
[31m-        """ return the keys of the variable dictionary"""[m
[31m-        return self.variableDict.keys()[m
[31m-    [m
     def save(self,index=None):[m
         """[m
         Saves the data dictionary to a NPY-file to a distinct location.[m
[1mdiff --git a/errors.py b/errors.py[m
[1mindex 834799c..fde714c 100644[m
[1m--- a/errors.py[m
[1m+++ b/errors.py[m
[36m@@ -7,7 +7,7 @@[m [mCreated on Fri Dec 14 11:44:59 2018[m
 """[m
 import inspect[m
 # =============================================================================[m
[31m-# [m
[32m+[m[32m#[m[41m  [m
 # =============================================================================[m
 class OperatorArgumentError(Exception):[m
     """Error if the operator has wrong arguments"""[m
[36m@@ -18,11 +18,19 @@[m [mclass OperatorArgumentError(Exception):[m
         msg = ''.join([msg1,msg2,msg3])[m
         super(OperatorArgumentError, self).__init__(msg)[m
 [m
[31m-def wrong_arguments(operator, correct_args):[m
[31m-    """ Throw an OperatorArgumentError if an operator does not have the correct arguments"""[m
[32m+[m[32mdef wrong_arguments(operator, expected_args):[m
[32m+[m[32m    """ Throws an OperatorArgumentError if an operator does not have the correct arguments/argument names[m[41m [m
[32m+[m[32m    needed for the module to be compatible with the from the user provided operators[m
[32m+[m[41m    [m
[32m+[m[32m    :param operator: the operator function that shell be tested for the expected arguments[m
[32m+[m[32m    :param expected_args: a list of strings of the expected arguments of the operator function[m
[32m+[m[41m    [m
[32m+[m[32m    :type operator: function[m
[32m+[m[32m    :type expected_args: list of strings[m
[32m+[m[32m    """[m
     operator_args = inspect.getargspec(operator).args[m
[31m-    if operator_args != correct_args:[m
[31m-        raise OperatorArgumentError(operator, operator_args, correct_args) [m
[32m+[m[32m    if operator_args != expected_args:[m
[32m+[m[32m        raise OperatorArgumentError(operator, operator_args, expected_args)[m[41m [m
 [m
 # =============================================================================[m
 #         [m
[36m@@ -34,7 +42,18 @@[m [mclass StateObjectKeyError(Exception):[m
         msg = "Used key: '%s' is not in list of initialised  keys: %s"%(test_key,correct_key_list)[m
         super(StateObjectKeyError, self).__init__(msg)[m
         [m
[31m-def check_keys(test_key_list, correct_key_list):[m
[32m+[m[32mdef check_keys(test_key_list, expected_key_list):[m
[32m+[m[32m    """ Throws an error if the key used for a state dictionary is not in the list of[m[41m [m
[32m+[m[32m    keys of the initialised stateObject (expected_key_list). Note, ref and tmp stateObject[m[41m [m
[32m+[m[32m    of the same category ("micro","macro","parameters") have the same 'expected_key_list'.[m
[32m+[m[41m    [m
[32m+[m[32m    :param test_key_list: list of keys that should be tested against the expected_key_list[m
[32m+[m[32m    :param expected_key_list: the expected_key_list[m
[32m+[m[41m    [m
[32m+[m[32m    :type test_key_list: list[m[41m [m
[32m+[m[32m    :type expected_key_list: list[m[41m [m
[32m+[m[32m    """[m
[32m+[m[41m    [m
     for test_key in test_key_list:[m
[31m-        if test_key not in correct_key_list:[m
[31m-            raise StateObjectKeyError(test_key, correct_key_list)[m
\ No newline at end of file[m
[32m+[m[32m        if test_key not in expected_key_list:[m
[32m+[m[32m            raise StateObjectKeyError(test_key, expected_key_list)[m
\ No newline at end of file[m
