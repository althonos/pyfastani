from . import test_ani, test_sketch

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_ani))
    suite.addTests(loader.loadTestsFromModule(test_sketch))
    return suite
