from . import test_mapper

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_mapper))
    return suite
