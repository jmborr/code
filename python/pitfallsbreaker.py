def pbhasattr(obj, name):
    """hasattr bypass class properties"""
    try:
        return getattr(obj, name, None) is not None
    except Exception, e:
        # Catches any error when name is a property being evaluated
        print "Error in property '{0}':".format(name), e
        return None
