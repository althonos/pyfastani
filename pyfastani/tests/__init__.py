from . import test_ani

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_ani))
    return suite
