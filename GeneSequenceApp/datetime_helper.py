def datetime_within_range(to_test, centre, margin):
        return (to_test - margin) < centre < (to_test + margin)
