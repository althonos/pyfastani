from . import test_ani, test_sketch, test_hasher

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_ani))
    suite.addTests(loader.loadTestsFromModule(test_sketch))
    suite.addTests(loader.loadTestsFromModule(test_hasher))
    return suite
