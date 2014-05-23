test_command_line:
	nosetests tests.test_montepython:Test01CommandLineInputBehaviour
test_setup:
	nosetests tests.test_montepython:Test02Setup
test_conf:
	nosetests tests.test_montepython:Test03NoDefaultConf
test_wrapper:
	nosetests tests.test_montepython:Test04CosmologicalCodeWrapper
test_data:
	nosetests tests.test_montepython:Test05DataModule
test_MH_IS:
	nosetests tests.test_montepython:Test06MetropolisHastingsImportanceSampling
test_CH:
	nosetests tests.test_montepython:Test07CosmoHammerBehaviour
test_NS:
	nosetests tests.test_montepython:Test08NestedSamplingBehaviour

short_tests: test_command_line test_setup test_conf test_CH
long_tests: test_wrapper test_data test_MH_IS test_NS

doctests:
	python montepython/data.py

tests: short_tests long_tests
all: doctests tests

