from rnaseqlyze.core.orm import Analysis, User

def create_anonymous_analysis(session):

    anonymous = session.query(User).get("anonymous")
    if not anonymous:
        anonymous = User("anonymous")
        session.commit()
    analysis = Analysis()
    analysis.owner = anonymous
    return analysis
