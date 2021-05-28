{% comment %}
As a maintainer, you don't need to edit this file.
If you notice that something doesn't work, please
open an issue: https://github.com/carpentries/styles/issues/new
{% endcomment %}

{% include manual_episode_order.html %}

<script>
  window.onload = function() {
    var lesson_episodes = [
    {% for lesson_episode in lesson_episodes %}
      {% if site.episode_order %}
        {% assign episode = site.episodes | where: "slug", lesson_episode | first %}
      {% else %}
        {% assign episode = lesson_episode %}
      {% endif %}
    "{{ episode.url}}"{% unless forloop.last %},{% endunless %}
    {% endfor %}
    ];
    var xmlHttp = [];  /* Required since we are going to query every episode. */
    for (i=0; i < lesson_episodes.length; i++) {
      xmlHttp[i] = new XMLHttpRequest();
      xmlHttp[i].episode = lesson_episodes[i];  /* To enable use this later. */
      xmlHttp[i].onreadystatechange = function() {
        if (this.readyState == 4 && this.status == 200) {
          var article_here = document.getElementById(this.episode);
          var parser = new DOMParser();
          var htmlDoc = parser.parseFromString(this.responseText,"text/html");
          var htmlDocArticle = htmlDoc.getElementsByTagName("article")[0];
          article_here.innerHTML = htmlDocArticle.innerHTML;
        }
      }
      var episode_url = "{{ relative_root_path }}" + lesson_episodes[i];
      xmlHttp[i].open("GET", episode_url);
      xmlHttp[i].send(null);
    }
  }
</script>

{% comment %} Create an anchor for every episode.  {% endcomment %}

{% for lesson_episode in lesson_episodes %}
  {% if site.episode_order %}
    {% assign episode = site.episodes | where: "slug", lesson_episode | first %}
  {% else %}
    {% assign episode = lesson_episode %}
  {% endif %}
  <article id="{{ episode.url }}"></article>
{% endfor %}
