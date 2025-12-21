import requests
import spotipy
from spotipy.oauth2 import SpotifyClientCredentials
import spotipy.util as util
import sys
import tidalapi

tidal_id = '68186599'
tidal_pwd = 'robercik2003_pl'
spotify_id = 'd5e63db4ea83471c'
tidal_username = 'robert.kow12321@gmail.com'
spotify_username = 'robert.kow12321@gmail.com'
SPOTIPY_CLIENT_ID = '9c93c9b33b7a41acacdf290722b40a0b'
SPOTIPY_CLIENT_SECRET = '5b89bf9be9514b88b2e8c1083d3dcee3'
SPOTIPY_REDIRECT_URI = 'http://127.0.0.1:8888/callback'
spotify_discover_weekly_id = 'your_spotify_discover_weekly_id'

# Options:
# connect_to_spotify()
# connect_to_tidal()
# delete_all_tidal_playlists()
# move_all_tidal_playlists_to_spotify()
# move_one_tidal_playlist_to_spotify(playlist_id)
# move_all_spotify_playlists_to_tidal()
# move_favourites_from_spotify_to_tidal()
# move_discover_weekly_from_spotify_to_tidal()

# Endpoints not in tidalapi
def get_tidal_create_playlist_url(tidal_id):
    return 'https://listen.tidal.com/v1/users/' + tidal_id + '/playlists'
def get_tidal_add_track_to_playlist_url(playlist_id):
    return 'https://listen.tidal.com/v1/playlists/' + playlist_id + '/items'
def get_tidal_find_track_url():
    return 'https://listen.tidal.com/v1/search/tracks'
def get_tidal_playlist(playlist_id):
    return 'https://listen.tidal.com/v1/playlists/' + playlist_id
def get_tidal_user_playlists():
    return 'https://listen.tidal.com/v1/users/' + tidal_id + '/playlists'

# Endpoints not in spotipy
def get_discover_weekly_playlist():
    return 'https://api.spotify.com/v1/users/spotify/playlists/'  + spotify_discover_weekly_id + '/tracks'

tidal_oldplaylists = []

def connect_to_spotify():
    scope = 'user-library-read'
    token = util.prompt_for_user_token(
        spotify_username,
        scope,
        client_id=SPOTIPY_CLIENT_ID,
        client_secret=SPOTIPY_CLIENT_SECRET,
        redirect_uri=SPOTIPY_REDIRECT_URI
    )

    if token:
        sp = spotipy.Spotify(auth=token)
    else:
        sys.exit()

    return sp, token


def connect_to_tidal():
    tidal_session = tidalapi.Session()
    try:
        tidal_session.login_oauth_simple()
    except requests.exceptions.HTTPError as e:
        print("Can't login to tidal for username=" + tidal_username + ", password=" + tidal_pwd)
        sys.exit()
    return tidal_session

sp, sp_token = connect_to_spotify()
tidal_session = connect_to_tidal()
        
def move_albums_from_spotify_to_tidal():
    user = tidal_session.user
    user_favourites = user.favorites

    def _get_spotify_favourites(offset):
        return sp.current_user_saved_albums(limit=20, offset=offset)

    albums = []
    offset = 0
    spotify_favourites = _get_spotify_favourites(offset)
    while spotify_favourites:
        for item in spotify_favourites['items']:
            # track = item['track']
            # track_id = _search_for_track_on_tidal(
            #     track['name'], track['artists'][0]['name']
            # )
            # if track_id > 0:
            #     user_favourites.add_track(track_id)
        offset = offset + 20
        spotify_favourites = _get_spotify_favourites(offset)


