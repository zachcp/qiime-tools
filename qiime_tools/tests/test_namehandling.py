import pytest
from qiime_tools.namehandling import processname


def test_namehandling():
    dotname = "test1.test2.test3.test4.test5.test6"
    assert(processname(name=dotname,func="dot") == "test1")
    assert(processname(name=dotname,func="dot2") == "test1.test2")
    assert(processname(name=dotname,func="dot3") == "test1.test2.test3")
    assert(processname(name=dotname,func="dot4") == "test1.test2.test3.test4")
    assert(processname(name=dotname,func="dot5") == "test1.test2.test3.test4.test5")

    underscorename = "test1_test2_test3_test4_test5_test6"
    assert(processname(name=underscorename,func="underscore") == "test1")
    assert(processname(name=underscorename,func="underscore2") == "test1.test2")
    assert(processname(name=underscorename,func="underscore3") == "test1.test2.test3")
    assert(processname(name=underscorename,func="underscore4") == "test1.test2.test3.test4")
    assert(processname(name=underscorename,func="underscore5") == "test1.test2.test3.test4.test5")


def test_namehandling_ErrorHandling():
    dotname = "test1.test2.test3.test4"
    with pytest.raises(ValueError) as error:
        processname(name=dotname, func="undeclaredfunc")
    assert 'You are attempting to process your sample names with an undefined function' in error.value.args[0]
