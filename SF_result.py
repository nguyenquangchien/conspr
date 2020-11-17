# SF_result.py
# Show the result from SF bay simulation

def getResult(par, dom):
    locations = [("Mateo Bridge", 75, 76), 
            ("East Oyster Pt.", 128, 76), 
            ("SW Oakland Airport", 113, 153),
            ("Site A", 132, 148),
            ("Site C", 116, 128),
            ]
    # Careful, index positions are reversed against X/Y
    if len(locations) < 6:
        for loc in locations:
            print loc[0], dom.h[loc[1], loc[2]], dom.u[loc[1], loc[2]], dom.v[loc[1], loc[2]]

